# tier2_hlae_fimbria_p53_vs_stic.R
# Same design as Tier 2 of neighborhood_metabolic_immune_analysis.R
# (HLA-E stratification within LDH+FASN double-positive index cells), but
# fit separately for TWO histology groups instead of STIC alone:
#   Group A: normalFimbria + p53signature combined ("Normal / p53sig")
#   Group B: STIC
# Incidental == Y only.
#
# Index cells: PanCK+, Incidental==Y, ROIcat %in% (normalFimbria,
#   p53signature, STIC), LDHp & FASNp (double-positive only)
# Groups within each roi_group: HLA-E- vs HLA-E+ (index cell's own status)
#
# GLMM per immune pop, per roi_group:
#   cbind(cnt, fails) ~ HLA_Ep_idx + (1 | slideName)
# Contrast: HLA-E+ vs HLA-E- (BH-FDR across 14 immune pops, WITHIN each
#   roi_group separately -- i.e. two independent correction families)
#
# NEW vs original Tier 2: per-slide balance filter -- within each roi_group
# independently, only slides with >= min_cells_per_slide cells in BOTH
# HLA-E- and HLA-E+ are kept before fitting (same threshold Tier 3 used).
#
# Outputs:
#   - one shared GLMM results CSV (both roi_groups)
#   - one shared summary-stats CSV (both roi_groups)
#   - two forest plots, one per roi_group (same visual style as
#     tier2_forest_plot.png: single dot + CI per population, colored only
#     by significance)

library(tidyverse)
library(arrow)
library(lme4)
library(emmeans)
library(scales)

# ── Settings ──────────────────────────────────────────────────────────────────
nbhd_parquet         <- "STIC_analysis/data_with_neighborhoods_50um_filt.parquet"
out_dir              <- "STIC_analysis/results/neighborhood_metabolic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_cells_per_slide  <- 5   # per-slide balance filter (HLA-E- and HLA-E+ each)

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

roi_levels <- c("Normal / p53sig", "STIC")

# ── Load & filter ─────────────────────────────────────────────────────────────
cat("Reading neighborhood parquet ...\n")
df_raw <- read_parquet(nbhd_parquet) %>%
  rename(HLA_Ep_idx = `HLA-Ep`)
cat(sprintf("  %d cells, %d slides\n", nrow(df_raw), n_distinct(df_raw$slideName)))

df <- df_raw %>%
  filter(
    PanCKp     == TRUE,
    Incidental == "Y",
    ROIcat %in% c("normalFimbria", "p53signature", "STIC"),
    LDHp       == TRUE,
    FASNp      == TRUE,
    !is.na(HLA_Ep_idx)
  ) %>%
  mutate(
    roi_group = case_when(
      ROIcat %in% c("normalFimbria", "p53signature") ~ "Normal / p53sig",
      ROIcat == "STIC"                                ~ "STIC"
    ),
    roi_group  = factor(roi_group, levels = roi_levels),
    HLA_Ep_idx = factor(HLA_Ep_idx, levels = c(FALSE, TRUE), labels = c("HLA-E-", "HLA-E+"))
  )

cat(sprintf("\nLDH+FASN index cells (both roi_groups, before slide filter): %d\n", nrow(df)))
cat("roi_group x HLA-E status:\n")
print(table(df$roi_group, df$HLA_Ep_idx))

# ── Per-slide balance filter (applied within each roi_group independently) ────
slide_counts <- df %>%
  count(slideName, roi_group, HLA_Ep_idx) %>%
  pivot_wider(names_from = HLA_Ep_idx, values_from = n, values_fill = 0L)

slides_ok <- slide_counts %>%
  filter(`HLA-E-` >= min_cells_per_slide, `HLA-E+` >= min_cells_per_slide) %>%
  distinct(slideName, roi_group)

cat(sprintf("\nSlides passing >=%d cells in both HLA-E- and HLA-E+:\n", min_cells_per_slide))
print(slides_ok %>% count(roi_group))

df <- df %>% semi_join(slides_ok, by = c("slideName", "roi_group"))
cat("\nCell counts after slide filter (used in GLMM):\n")
print(table(df$roi_group, df$HLA_Ep_idx))
cat("Slides used per roi_group:\n")
print(df %>% group_by(roi_group) %>% summarise(n_slides = n_distinct(slideName)))

# ── Summary statistics table ───────────────────────────────────────────────────
# Per roi_group x HLA-E status: n cells, n slides, and median/mean % of each
# immune population within the 50um neighborhood (mirrors Tier 2 boxplots).
summary_counts <- df %>%
  group_by(roi_group, HLA_Ep_idx) %>%
  summarise(n_cells = n(), n_slides = n_distinct(slideName), .groups = "drop")

summary_pct <- df %>%
  pivot_longer(cols = all_of(immune_pct_cols), names_to = "immune_col", values_to = "pct") %>%
  mutate(immune_pop = sub("_pct$", "", immune_col)) %>%
  group_by(roi_group, HLA_Ep_idx, immune_pop) %>%
  summarise(
    median_pct = median(pct, na.rm = TRUE),
    mean_pct   = mean(pct, na.rm = TRUE),
    .groups    = "drop"
  )

summary_stats <- summary_pct %>%
  left_join(summary_counts, by = c("roi_group", "HLA_Ep_idx")) %>%
  mutate(immune_pop = factor(immune_pop, levels = immune_pops)) %>%
  arrange(roi_group, immune_pop, HLA_Ep_idx) %>%
  mutate(across(c(median_pct, mean_pct), ~ round(.x, 2)))

write_csv(summary_stats, file.path(out_dir, "tier2_hlae_fimbria_p53_vs_stic_summary_stats.csv"))
cat("\n=== Summary stats (first 20 rows) ===\n")
print(summary_stats, n = 20)

# ── GLMM helper: one immune population, one roi_group ──────────────────────────
run_glmm <- function(cnt_col, data) {
  pop        <- sub("_count$", "", cnt_col)
  data$cnt   <- data[[cnt_col]]
  data$fails <- data$total_neighbors - data$cnt
  data        <- filter(data, fails >= 0, !is.na(cnt))

  empty <- tibble(contrast = NA_character_, odds.ratio = NA_real_,
                  asymp.LCL = NA_real_, asymp.UCL = NA_real_,
                  p.value = NA_real_, immune_pop = pop, note = "skipped")

  if (nrow(data) < 30 || n_distinct(data$slideName) < 2 ||
      n_distinct(data$HLA_Ep_idx) < 2)
    return(mutate(empty, note = "insufficient_data"))

  form <- as.formula(paste0(
    "cbind(`", cnt_col, "`, total_neighbors - `", cnt_col, "`) ~ ",
    "HLA_Ep_idx + (1 | slideName)"
  ))
  m <- tryCatch(
    glmer(form, family = binomial, data = data,
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
    error = function(e) NULL
  )
  if (is.null(m)) return(mutate(empty, note = "model_failed"))

  em <- tryCatch(emmeans(m, ~ HLA_Ep_idx), error = function(e) NULL)
  if (is.null(em)) return(mutate(empty, note = "emmeans_failed"))

  contrast(em, "trt.vs.ctrl", ref = 1, type = "response", infer = c(TRUE, TRUE)) %>%
    as_tibble() %>%
    mutate(contrast = "HLA-E+ vs HLA-E-", immune_pop = pop, note = "ok")
}

run_glmm_for_group <- function(group_name, data) {
  cat(sprintf("\nFitting GLMMs for roi_group = %s (14 pops) ...\n", group_name))
  grp_data <- data %>% filter(roi_group == group_name)
  map_dfr(immune_cnt_cols, run_glmm, data = grp_data) %>%
    mutate(roi_group = group_name)
}

res_all <- map_dfr(roi_levels, run_glmm_for_group, data = df) %>%
  filter(!is.na(odds.ratio)) %>%
  mutate(immune_pop = factor(immune_pop, levels = immune_pops)) %>%
  group_by(roi_group) %>%
  mutate(
    p_adj     = p.adjust(p.value, method = "BH"),
    sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    significant = p_adj < 0.05
  ) %>%
  ungroup() %>%
  arrange(roi_group, p_adj)

cat("\n=== Results (HLA-E+ vs HLA-E-, both roi_groups, BH-FDR WITHIN each group) ===\n")
res_all %>%
  select(roi_group, immune_pop, odds.ratio, asymp.LCL, asymp.UCL, p.value, p_adj, sig_label) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print(n = Inf)

write_csv(res_all, file.path(out_dir, "tier2_hlae_fimbria_p53_vs_stic_glmm_results.csv"))

# ── Forest plots: one file per roi_group ────────────────────────────────────────
# Same visual convention as tier2_forest_plot.png: one dot + CI per
# population, colored only by significance.

make_forest <- function(data, group_name) {
  plot_dat <- data %>%
    filter(roi_group == group_name, !is.na(odds.ratio)) %>%
    mutate(immune_pop = factor(as.character(immune_pop), levels = rev(immune_pops)))

  ggplot(plot_dat, aes(x = odds.ratio, y = immune_pop,
                       color = significant, alpha = significant)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_point(size = 2.8) +
    geom_errorbar(aes(xmin = asymp.LCL, xmax = asymp.UCL),
                  width = 0.45, linewidth = 0.6, orientation = "y") +
    geom_text(
      data = function(d) filter(d, significant),
      aes(x = asymp.UCL, label = sig_label),
      hjust = -0.25, size = 3.2, show.legend = FALSE
    ) +
    scale_x_log10(name = "Odds ratio (log scale): HLA-E+ vs HLA-E-") +
    scale_color_manual(
      values = c("TRUE" = "#D73027", "FALSE" = "grey55"),
      labels = c("TRUE" = "BH adj. p < 0.05", "FALSE" = "ns"),
      name   = NULL
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), guide = "none") +
    labs(
      y = NULL,
      title = str_wrap(paste0(group_name, ": Immune neighborhood — HLA-E+ vs HLA-E- within LDH+FASN cells"), width = 62),
      subtitle = paste("PanCK+  |  Incidental==Y  |  LDHp & FASNp  |",
                       sprintf(">=%d cells/slide/group  |  BH-FDR", min_cells_per_slide))
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold"),
      axis.text.y   = element_text(size = 9),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      legend.position = "right",
      plot.margin = margin(10, 20, 10, 10)
    )
}

for (grp in roi_levels) {
  grp_slug <- if (grp == "Normal / p53sig") "fimbria_p53" else "stic"
  p <- make_forest(res_all, grp)

  ggsave(file.path(out_dir, paste0("tier2_forest_", grp_slug, ".pdf")), p, width = 10, height = 7)
  ggsave(file.path(out_dir, paste0("tier2_forest_", grp_slug, ".png")), p, width = 10, height = 7, dpi = 200)
  cat(sprintf("  Saved forest plot for: %s\n", grp))
}

message("Done. Outputs saved to: ", out_dir)
