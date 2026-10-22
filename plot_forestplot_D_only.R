library(tidyverse)

roi_D <- c("normalFimbria + p53sig", "STIL", "STIC")

res_D <- read_csv("Ki67_LDH_FISHgroup_comparison/FISH_Ki67_LDH_D_GLMM_results.csv",
                  show_col_types = FALSE)

# OR range: 0.507 – 2.35; CI range: 0.105 – 4.34 (STIL CCNE1 only, both markers)
# Breaks span the full CI range; no clipping applied
or_breaks_D <- c(0.1, 0.2, 0.5, 1.0, 2.0, 4.0)

p <- res_D |>
  mutate(
    FISH_group   = factor(FISH_group, levels = c("MYC only", "CCNE1 only", "Both")),
    ROIcat       = factor(ROIcat, levels = roi_D),
    outcome      = recode(outcome,
      Ki67 = "Ki67 (Proliferation)",
      LDH  = "LDH (Metabolism)"
    ),
    significance = replace(significance, is.na(significance), "")
  ) |>
  ggplot(aes(x = OR, y = FISH_group, color = FISH_group)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_errorbar(aes(xmin = OR_ci_lo, xmax = OR_ci_hi),
                width = 0.25, linewidth = 0.8, orientation = "y") +
  geom_point(size = 3.5) +
  geom_text(aes(label = significance, x = OR_ci_hi),
            hjust = -0.3, size = 4.5, color = "black", fontface = "bold") +
  scale_color_manual(
    values = c("MYC only" = "#4393c3", "CCNE1 only" = "#d6604d", "Both" = "#b2182b"),
    name   = "FISH group"
  ) +
  scale_x_log10(
    breaks = or_breaks_D,
    labels = or_breaks_D,
    expand = expansion(mult = c(0.05, 0.20))
  ) +
  facet_grid(ROIcat ~ outcome) +
  labs(
    title    = "Odds of Ki67+ and LDH+ vs Diploid, by ROI Type",
    subtitle = "Incidental cases — E-cadherin+ cells — GLMM (1 | slideName)\nnormalFimbria + p53sig  |  STIL  |  STIC",
    x        = "Odds ratio vs Diploid  (log scale, ± 95% CI)",
    y        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10),
    strip.text       = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    axis.text.x      = element_text(size = 10),
    axis.text.y      = element_text(size = 10),
    legend.position  = "bottom",
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey70", fill = NA)
  )

ggsave("FISH_Ki67_LDH_OR_forestplot_D.pdf", p, width = 13, height = 9)
ggsave("FISH_Ki67_LDH_OR_forestplot_D.png", p, width = 13, height = 9, dpi = 200)
