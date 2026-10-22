# neighborhood_metabolic_immune_analysis.R
# Reads pre-computed 50um neighborhood parquet (from build_nbhd_50um.R)
#
# TIER 1: Metabolic neighborhood comparison
#   Index cells: PanCK+, Incidental==Y, normalFimbria/p53sig/STIL/STIC
#   4 groups: Neither / LDH only / FASN only / LDH+FASN
#   GLMM: cbind(cnt, fails) ~ LDHp_f * FASNp_f + (1|slideName), per ROIcat group
#   3 contrasts vs Neither (BH-FDR across 14 immune pops per contrast x ROI group)


library(tidyverse)
library(arrow)
library(lme4)
library(emmeans)
library(broom.mixed)
library(scales)
library(forcats)

# ── Settings ──────────────────────────────────────────────────────────────────
nbhd_parquet    <- "STIC_analysis/data_with_neighborhoods_50um_filt.parquet"
out_dir         <- "STIC_analysis/results/neighborhood_metabolic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_idx_per_roi <- 500   # index cells subsampled per (slide x ROI) for GLMM speed
set.seed(42)

immune_pops <- c(
  "CD8ap", "CD4p", "NKG2ApCD8anCD4n", "CD11cp", "CD68p",
  "CD8apLAG3pTIM3p", "CD8apTbetpNKG2ApCD103pLAG3pTIM3p",
  "CD8apTbetpNKG2ApCD103nLAG3pTIM3p", "CD4pTIM3pLAG3p",
  "CD11cpCD68pCD103nTIM3p", "CD11cpCD68nCD103nTIM3p",
  "CD8apTbetpNKG2ApCD103pPD1pKi67pLAG3nTIM3n",
  "CD8apTbetpNKG2ApCD103nPD1pKi67pLAG3nTIM3n", "CD68pCD11cp"
)
immune_cnt_cols <- paste0(immune_pops, "_count")
immune_pct_cols <- paste0(immune_pops, "_pct")

roi_levels <- c("Normal / p53sig", "STIL", "STIC")

group_pal <- c(
  "Neither"   = "#808080",
  "LDH only"  = "#E69F00",
  "FASN only" = "#56B4E9",
  "LDH+FASN"  = "#D55E00"
)
contrast_pal <- c(
  "LDH only vs Neither"  = "#E69F00",
  "FASN only vs Neither" = "#56B4E9",
  "LDH+FASN vs Neither"  = "#D55E00"
)
hlae_pal <- c("HLA-E-" = "#4575B4", "HLA-E+" = "#D73027")

base_theme <- function() {
  theme_bw(base_size = 10) +
    theme(
      strip.text       = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

# ── Load parquet ──────────────────────────────────────────────────────────────
cat("Reading neighborhood parquet ...\n")
df_raw <- read_parquet(nbhd_parquet) %>%
  rename(HLA_Ep_idx = `HLA-Ep`)   # index cell own HLA-E status; avoid hyphen in formulas

cat(sprintf("  %d cells, %d slides\n", nrow(df_raw), n_distinct(df_raw$slideName)))
cat("  ROIcat distribution:\n")
print(table(df_raw$ROIcat, useNA = "ifany"))

# =============================================================================
# TIER 1: Metabolic neighborhood comparison
# =============================================================================
cat("\n====== TIER 1: Metabolic neighborhood comparison ======\n")

df_t1 <- df_raw %>%
  filter(
    PanCKp     == TRUE,
    Incidental == "Y",
    ROIcat %in% c("normalFimbria", "p53signature", "STIC", "STIL")
  ) %>%
  mutate(
    roi_group = case_when(
      ROIcat %in% c("normalFimbria", "p53signature") ~ "Normal / p53sig",
      ROIcat == "STIL"                               ~ "STIL",
      ROIcat == "STIC"                               ~ "STIC"
    ),
    roi_group      = factor(roi_group, levels = roi_levels),
    LDH_FASN_group = factor(LDH_FASN_group,
                            levels = c("Neither","LDH only","FASN only","LDH+FASN")),
    LDHp_f  = factor(LDHp,  levels = c(FALSE, TRUE), labels = c("LDHn", "LDHp")),
    FASNp_f = factor(FASNp, levels = c(FALSE, TRUE), labels = c("FASNn", "FASNp"))
  )

cat(sprintf("  Cells: %d\n", nrow(df_t1)))
cat("  roi_group x LDH_FASN_group:\n")
print(table(df_t1$roi_group, df_t1$LDH_FASN_group))

# Subsample index cells per (slide x ROI x roi_group) for GLMM speed
df_t1_sub <- df_t1 %>%
  group_by(slideName, ROIname, roi_group) %>%
  slice_sample(n = max_idx_per_roi, replace = FALSE) %>%
  ungroup()
cat(sprintf("  Subsampled for GLMM: %d cells\n", nrow(df_t1_sub)))

# ── GLMM helper: one immune pop, all ROI groups ───────────────────────────────
run_t1_glmm <- function(cnt_col, data) {
  pop <- sub("_count$", "", cnt_col)
  data %>%
    group_by(roi_group) %>%
    group_modify(function(grp, key) {
      grp$cnt   <- grp[[cnt_col]]
      grp$fails <- grp$total_neighbors - grp$cnt
      grp        <- filter(grp, fails >= 0, !is.na(cnt))

      empty <- tibble(contrast = NA_character_, odds.ratio = NA_real_,
                      asymp.LCL = NA_real_, asymp.UCL = NA_real_,
                      p.value = NA_real_, immune_pop = pop, note = "skipped")

      if (nrow(grp) < 30 || n_distinct(grp$slideName) < 2) return(empty)

      form <- as.formula(paste0(
        "cbind(`", cnt_col, "`, total_neighbors - `", cnt_col, "`) ~ ",
        "LDHp_f * FASNp_f + (1 | slideName)"
      ))
      m <- tryCatch(
        glmer(form, family = binomial, data = grp,
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl   = list(maxfun = 2e5))),
        error = function(e) NULL
      )
      if (is.null(m)) return(mutate(empty, note = "model_failed"))

      em <- tryCatch(emmeans(m, ~ LDHp_f * FASNp_f), error = function(e) NULL)
      if (is.null(em)) return(mutate(empty, note = "emmeans_failed"))

      contrast(em, "trt.vs.ctrl", ref = 1, type = "response",
               infer = c(TRUE, TRUE)) %>%
        as_tibble() %>%
        mutate(
          contrast = case_when(
            grepl("LDHn FASNp", contrast) ~ "FASN only vs Neither",
            grepl("LDHp FASNn", contrast) ~ "LDH only vs Neither",
            grepl("LDHp FASNp", contrast) ~ "LDH+FASN vs Neither",
            TRUE ~ contrast
          ),
          immune_pop = pop,
          note = "ok"
        )
    }) %>%
    ungroup()
}

cat("  Fitting Tier 1 GLMMs (14 pops x 3 ROI groups) ...\n")
t1_res <- map_dfr(immune_cnt_cols, run_t1_glmm, data = df_t1_sub)

# BH-FDR within each (roi_group x contrast) across 14 immune pops
t1_res <- t1_res %>%
  filter(!is.na(odds.ratio)) %>%
  group_by(roi_group, contrast) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    sig_label   = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                            p_adj < 0.05  ~ "*",   TRUE ~ "ns"),
    significant = p_adj < 0.05,
    immune_pop  = factor(immune_pop, levels = rev(immune_pops)),
    contrast    = factor(contrast, levels = c(
      "LDH only vs Neither", "FASN only vs Neither", "LDH+FASN vs Neither"
    ))
  )

write.csv(t1_res, file.path(out_dir, "tier1_glmm_results.csv"), row.names = FALSE)

cat("\n  Significant Tier 1 results (BH adj p < 0.05):\n")
t1_res %>%
  filter(significant) %>%
  arrange(roi_group, contrast, p_adj) %>%
  select(roi_group, immune_pop, contrast, odds.ratio, p.value, p_adj, sig_label) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print(n = Inf)

# ── Tier 1 forest plots (split: basic immune populations vs cell states) ──────

basic_pops <- c("CD8ap", "CD4p", "NKG2ApCD8anCD4n", "CD11cp", "CD68p", "CD68pCD11cp")
state_pops <- c(
  "CD8apLAG3pTIM3p", "CD8apTbetpNKG2ApCD103pLAG3pTIM3p",
  "CD8apTbetpNKG2ApCD103nLAG3pTIM3p", "CD4pTIM3pLAG3p",
  "CD11cpCD68pCD103nTIM3p", "CD11cpCD68nCD103nTIM3p",
  "CD8apTbetpNKG2ApCD103pPD1pKi67pLAG3nTIM3n",
  "CD8apTbetpNKG2ApCD103nPD1pKi67pLAG3nTIM3n"
)

make_t1_forest <- function(data, pop_subset, title_suffix) {
  plot_dat <- data %>%
    filter(!is.na(odds.ratio),
           as.character(immune_pop) %in% pop_subset) %>%
    mutate(immune_pop = factor(as.character(immune_pop),
                               levels = rev(pop_subset)))
  sig_dat <- filter(plot_dat, significant)

  ggplot(plot_dat, aes(x = odds.ratio, y = immune_pop,
                       color = contrast, alpha = significant)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_point(position = position_dodge(0.65), size = 2.3) +
    geom_errorbar(
      aes(xmin = asymp.LCL, xmax = asymp.UCL),
      position  = position_dodge(0.65),
      width     = 0.35, linewidth = 0.5,
      orientation = "y"
    ) +
    geom_text(
      data = sig_dat,
      aes(x = asymp.UCL, label = sig_label),
      position = position_dodge(0.65),
      hjust = -0.2, size = 2.5, show.legend = FALSE
    ) +
    scale_x_log10(name = "Odds ratio (log scale) vs Neither") +
    scale_color_manual(values = contrast_pal, name = "Contrast") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.22), guide = "none") +
    facet_wrap(~ roi_group, ncol = 3) +
    labs(
      y = NULL,
      title    = paste("Tier 1:", title_suffix),
      subtitle = "PanCK+  |  Incidental==Y  |  OR vs Neither (LDH-FASN-)  |  BH-FDR across 14 pops"
    ) +
    base_theme()
}

p_t1_forest_basic <- make_t1_forest(
  t1_res, basic_pops,
  "Immune neighborhood — basic immune populations"
)
ggsave(file.path(out_dir, "tier1_forest_plot_basic.pdf"),
       p_t1_forest_basic, width = 13, height = 5)
ggsave(file.path(out_dir, "tier1_forest_plot_basic.png"),
       p_t1_forest_basic, width = 13, height = 5, dpi = 180)

p_t1_forest_states <- make_t1_forest(
  t1_res, state_pops,
  "Immune neighborhood — immune cell states"
)
ggsave(file.path(out_dir, "tier1_forest_plot_states.pdf"),
       p_t1_forest_states, width = 13, height = 7)
ggsave(file.path(out_dir, "tier1_forest_plot_states.png"),
       p_t1_forest_states, width = 13, height = 7, dpi = 180)

# Heatmap: median % neighbors, all ROI groups
df_t1_heat <- df_t1 %>%
  pivot_longer(cols = all_of(immune_pct_cols),
               names_to = "immune_col", values_to = "pct") %>%
  mutate(immune_pop = factor(sub("_pct$", "", immune_col), levels = rev(immune_pops))) %>%
  group_by(roi_group, LDH_FASN_group, immune_pop) %>%
  summarise(median_pct = median(pct, na.rm = TRUE), .groups = "drop")

p_t1_heat <- ggplot(
  df_t1_heat,
  aes(x = LDH_FASN_group, y = immune_pop, fill = median_pct)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", name = "Median %\nin 50um nbhd") +
  facet_wrap(~ roi_group) +
  labs(x = NULL, y = NULL,
       title = "Tier 1: Median % immune neighbors by metabolic group",
       subtitle = "PanCK+  |  Incidental==Y") +
  base_theme() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave(file.path(out_dir, "tier1_heatmap.pdf"),
       p_t1_heat, width = 13, height = 7)
ggsave(file.path(out_dir, "tier1_heatmap.png"),
       p_t1_heat, width = 13, height = 7, dpi = 180)

# Boxplots — STIC only for focus
p_t1_box <- df_t1 %>%
  filter(roi_group == "STIC") %>%
  pivot_longer(cols = all_of(immune_pct_cols),
               names_to = "immune_col", values_to = "pct") %>%
  mutate(immune_pop = factor(sub("_pct$", "", immune_col), levels = immune_pops)) %>%
  ggplot(aes(x = LDH_FASN_group, y = pct, fill = LDH_FASN_group)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3,
               alpha = 0.8, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 0.35, alpha = 0.25) +
  scale_fill_manual(values = group_pal, guide = "none") +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  facet_wrap(~ immune_pop, scales = "free_y", ncol = 4) +
  labs(x = NULL, y = "% of 50um neighbors",
       title = "Tier 1: Immune neighborhood — STIC (Incidental==Y)",
       subtitle = "Each point = one PanCK+ index cell") +
  base_theme() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7))

ggsave(file.path(out_dir, "tier1_boxplots_STIC.pdf"),
       p_t1_box, width = 14, height = 11)
ggsave(file.path(out_dir, "tier1_boxplots_STIC.png"),
       p_t1_box, width = 14, height = 11, dpi = 180)

message("Tier 1 outputs saved.")

