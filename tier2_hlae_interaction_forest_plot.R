# tier2_hlae_interaction_forest_plot.R
# Forest plot of the HLA-E x roi_group interaction test results produced by
# tier2_hlae_interaction_test.R (reads the saved CSV, no refitting).
#
# Shows the ratio of odds ratios (ROR = OR[STIC] / OR[Normal/p53sig]) per
# immune population, on a log scale centered at ROR = 1 (no group
# difference). Colored by direction of the significant difference (rather
# than just significant/ns, since ROR can move in either direction).

library(tidyverse)
library(scales)

out_dir <- "STIC_analysis/results/neighborhood_metabolic"

immune_pops <- c(
  "CD8ap", "CD4p", "NKG2ApCD8anCD4n", "CD11cp", "CD68p",
  "CD8apLAG3pTIM3p", "CD8apTbetpNKG2ApCD103pLAG3pTIM3p",
  "CD8apTbetpNKG2ApCD103nLAG3pTIM3p", "CD4pTIM3pLAG3p",
  "CD11cpCD68pCD103nTIM3p", "CD11cpCD68nCD103nTIM3p",
  "CD8apTbetpNKG2ApCD103pPD1pKi67pLAG3nTIM3n",
  "CD8apTbetpNKG2ApCD103nPD1pKi67pLAG3nTIM3n", "CD68pCD11cp"
)

ror_res <- read_csv(file.path(out_dir, "tier2_hlae_interaction_ROR_results.csv"),
                     show_col_types = FALSE) %>%
  mutate(
    immune_pop = factor(immune_pop, levels = rev(immune_pops)),
    direction  = factor(
      direction,
      levels = c("significantly stronger in STIC",
                 "significantly stronger in Normal/p53sig",
                 "no significant difference")
    )
  )

dir_pal <- c(
  "significantly stronger in STIC"          = "#D73027",
  "significantly stronger in Normal/p53sig" = "#4575B4",
  "no significant difference"               = "grey55"
)

p <- ggplot(ror_res, aes(x = ROR, y = immune_pop, color = direction, alpha = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(xmin = ROR.LCL, xmax = ROR.UCL),
                width = 0.45, linewidth = 0.6, orientation = "y") +
  geom_text(
    data = function(d) filter(d, significant),
    aes(x = ROR.UCL, label = sig_label),
    hjust = -0.25, size = 3.2, show.legend = FALSE
  ) +
  scale_x_log10(
    name   = "Ratio of odds ratios (ROR): HLA-E effect in STIC / HLA-E effect in Normal-p53sig",
    breaks = c(0.5, 0.6, 0.7, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)
  ) +
  scale_color_manual(values = dir_pal, name = NULL) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), guide = "none") +
  labs(
    y = NULL,
    title = str_wrap("Interaction test: does the HLA-E+ vs HLA-E- neighborhood effect differ between Normal/p53sig and STIC?", width = 70),
    subtitle = paste(
      "PanCK+  |  Incidental==Y  |  LDHp & FASNp  |  >=5 cells/slide/group  |  BH-FDR\n",
      "cbind(count, fails) ~ HLA-E status x histology group + (1|slideName); ROR > 1 = stronger effect in STIC"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 9),
    plot.subtitle   = element_text(size = 7.5, color = "grey40"),
    legend.position = "right",
    plot.margin     = margin(10, 20, 10, 10)
  )

ggsave(file.path(out_dir, "tier2_hlae_interaction_forest_plot.pdf"), p, width = 11, height = 7)
ggsave(file.path(out_dir, "tier2_hlae_interaction_forest_plot.png"), p, width = 11, height = 7, dpi = 200)

message("Saved: tier2_hlae_interaction_forest_plot.pdf / .png")
