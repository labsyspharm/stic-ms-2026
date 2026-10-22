# tier3_triple_vs_neither_STIC.R
# Compares 50um immune neighborhood of Triple+ (LDH+FASN+Ki67+) vs
# Neither (LDH-FASN-Ki67-) PanCK+ tumor cells in STIC (Incidental==Y).
# Uses pre-computed neighborhoods from data_with_neighborhoods_50um_filt.parquet.

library(tidyverse)
library(arrow)
library(lme4)
library(emmeans)
library(scales)

input_parquet <- "STIC_analysis/data_with_neighborhoods_50um_filt.parquet"
out_dir       <- "STIC_analysis/results/neighborhood_metabolic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

immune_pops <- c(
  "CD8ap", "CD4p", "NKG2ApCD8anCD4n", "CD11cp", "CD68p",
  "CD8apLAG3pTIM3p", "CD8apTbetpNKG2ApCD103pLAG3pTIM3p",
  "CD8apTbetpNKG2ApCD103nLAG3pTIM3p", "CD4pTIM3pLAG3p",
  "CD11cpCD68pCD103nTIM3p", "CD11cpCD68nCD103nTIM3p",
  "CD8apTbetpNKG2ApCD103pPD1pKi67pLAG3nTIM3n",
  "CD8apTbetpNKG2ApCD103nPD1pKi67pLAG3nTIM3n", "CD68pCD11cp"
)
cnt_cols <- paste0(immune_pops, "_count")
pct_cols <- paste0(immune_pops, "_pct")

# ── Load & filter ─────────────────────────────────────────────────────────────
cat("Reading neighborhood parquet ...\n")
dat <- read_parquet(input_parquet)
cat(sprintf("  %d cells, %d columns\n", nrow(dat), ncol(dat)))

t3 <- dat |>
  filter(ROIcat == "STIC", PanCKp == TRUE, Incidental == "Y") |>
  mutate(
    metabolic_group = case_when(
       LDHp &  FASNp &  Ki67p ~ "Triple",
      !LDHp & !FASNp & !Ki67p ~ "Neither",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(metabolic_group), total_neighbors > 0) |>
  mutate(metabolic_group = factor(metabolic_group, levels = c("Neither", "Triple")))

cat("\n=== Tier 3: Triple vs Neither in STIC (PanCK+, Incidental==Y) ===\n")
cat("Group counts:\n")
print(table(t3$metabolic_group))

# Keep only slides with >= 5 cells in BOTH groups
slide_counts <- t3 |>
  count(slideName, metabolic_group) |>
  pivot_wider(names_from = metabolic_group, values_from = n, values_fill = 0L)
slides_ok <- slide_counts |> filter(Neither >= 5, Triple >= 5)
cat(sprintf("\nSlides with >=5 cells in both groups: %d\n", nrow(slides_ok)))

t3 <- t3 |> filter(slideName %in% slides_ok$slideName)
cat("After slide filter:\n")
print(table(t3$metabolic_group))
cat("Per slide:\n")
print(table(t3$slideName, t3$metabolic_group))

if (nrow(slides_ok) < 3) stop("Too few informative slides for GLMM.")

# ── GLMM per immune population ────────────────────────────────────────────────
cat("\nFitting Tier 3 GLMMs (14 pops, Triple vs Neither, STIC) ...\n")

run_t3_glmm <- function(cnt_col) {
  pop <- sub("_count$", "", cnt_col)
  tryCatch({
    grp <- t3 |>
      select(slideName, metabolic_group, total_neighbors,
             cnt = all_of(cnt_col)) |>
      filter(!is.na(cnt)) |>
      mutate(fails = total_neighbors - cnt)

    m <- glmer(
      cbind(cnt, fails) ~ metabolic_group + (1 | slideName),
      family  = binomial,
      data    = grp,
      control = glmerControl(optimizer = "bobyqa",
                             optCtrl   = list(maxfun = 2e5))
    )
    em  <- emmeans(m, ~ metabolic_group)
    res <- contrast(em, "trt.vs.ctrl", ref = 1, type = "response",
                    infer = c(TRUE, TRUE)) |>
      as_tibble() |>
      mutate(immune_pop = pop, contrast = "Triple vs Neither")
    res
  }, error = function(e) {
    warning("GLMM failed for ", pop, ": ", conditionMessage(e))
    NULL
  })
}

t3_results <- map(cnt_cols, run_t3_glmm) |>
  compact() |>
  bind_rows() |>
  mutate(
    immune_pop = factor(immune_pop, levels = immune_pops),
    p_adj      = p.adjust(p.value, method = "BH"),
    sig_label  = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  ) |>
  arrange(p_adj)

cat("\nTier 3 results (Triple vs Neither, STIC, BH-FDR):\n")
print(t3_results |>
  select(immune_pop, odds.ratio, asymp.LCL, asymp.UCL, p.value, p_adj, sig_label))

write_csv(t3_results, file.path(out_dir, "tier3_glmm_results.csv"))

# ── Forest plot ────────────────────────────────────────────────────────────────
forest_dat <- t3_results |>
  mutate(
    immune_pop = fct_reorder(immune_pop, odds.ratio),
    sig        = sig_label != "ns"
  )

p_forest <- ggplot(forest_dat,
    aes(x = odds.ratio, y = immune_pop, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = asymp.LCL, xmax = asymp.UCL),
                orientation = "y", width = 0.4, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(label = sig_label, x = asymp.UCL),
            hjust = -0.3, size = 3.5, color = "black") +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#8E44AD"),
                     guide = "none") +
  scale_x_log10(labels = label_number(accuracy = 0.1)) +
  labs(
    x        = "Odds Ratio (Triple vs Neither, log scale)",
    y        = NULL,
    title    = "Immune neighborhood: Triple+ vs Neither (STIC)",
    subtitle = "LDH+FASN+Ki67+ vs LDH-FASN-Ki67-  |  PanCK+, Incidental  |  50 um radius  |  GLMM binomial + BH-FDR"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold"),
    axis.text.y   = element_text(size = 9),
    plot.subtitle = element_text(size = 8, color = "grey40")
  )

ggsave(file.path(out_dir, "tier3_forest_plot.pdf"), p_forest, width = 9, height = 6)
ggsave(file.path(out_dir, "tier3_forest_plot.png"), p_forest, width = 9, height = 6, dpi = 200)
cat("Forest plot saved.\n")

