## Per-patient Immune Surveillance Index, stratified by the JOINT IFN x Metabolic
## tertile state of the STIC epithelial index cell (4 groups: Low_Low, Low_High,
## High_Low, High_High -- Mid tertile excluded on both axes). Extends
## STIC_immune_surveillance_score.R, which only stratified by one axis at a time.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

data_dir <- "c:/Users/tak415/OneDrive - Harvard University/Documents/CyCIF_Geomx_analysis_LSP_HMS/CYCIF_ORION_atlas_2.0/Xenium_analysis"
out_dir  <- data_dir

ifn <- readRDS(file.path(data_dir, "STIC_epithelial_neighborhood_readouts.rds"))
met <- readRDS(file.path(data_dir, "STIC_epithelial_metabolic_neighborhood_readouts.rds"))

shared_cols <- c("frac_CD8","frac_CD4","frac_NK","frac_MacDC","frac_HAVCR2_of_MacDC",
                  "CD8_exhaustion","CD8_activation","CD4_exhaustion","CD4_activation","NK_activation")

unified <- ifn %>%
  select(cell, sample_id, IFN_tertile, all_of(shared_cols)) %>%
  inner_join(met %>% select(cell, sample_id, Metab_tertile), by = c("cell","sample_id"))

## 1. Joint group per cell (Low/High only on both axes) ----------------------
j <- unified %>%
  filter(IFN_tertile %in% c("Low","High"), Metab_tertile %in% c("Low","High")) %>%
  mutate(joint_group = paste0(IFN_tertile, "IFN_", Metab_tertile, "Metab"))
j$joint_group <- factor(j$joint_group,
                         levels = c("LowIFN_LowMetab", "LowIFN_HighMetab", "HighIFN_LowMetab", "HighIFN_HighMetab"))

density_cols <- c("frac_CD8","frac_CD4","frac_NK","frac_MacDC")

## 2. Per-patient x joint-group means, with cell-count feasibility flag ------
pm <- j %>%
  group_by(sample_id, joint_group) %>%
  summarise(n_cells = n(),
            across(all_of(c(density_cols, "frac_HAVCR2_of_MacDC",
                             "CD8_activation","CD8_exhaustion","CD4_activation","CD4_exhaustion","NK_activation")),
                   ~mean(.x, na.rm = TRUE)),
            .groups = "drop")
cat("Cell counts per patient x joint-group (flagging n<5):\n")
print(pm %>% select(sample_id, joint_group, n_cells) %>% mutate(low_n_flag = n_cells < 5) %>%
        pivot_wider(id_cols = sample_id, names_from = joint_group, values_from = n_cells))

## 3. Z-score across all 56 rows jointly, build Density_z / Function_z / Surveillance
z <- function(x) as.numeric(scale(x))
pm <- pm %>% mutate(across(all_of(c(density_cols, "CD8_activation","CD8_exhaustion",
                                     "CD4_activation","CD4_exhaustion","NK_activation",
                                     "frac_HAVCR2_of_MacDC")),
                            z, .names = "z_{.col}"))
pm <- pm %>% mutate(
  Density_z  = rowMeans(across(all_of(paste0("z_", density_cols))), na.rm = TRUE),
  CD8_net_z  = z_CD8_activation - z_CD8_exhaustion,
  CD4_net_z  = z_CD4_activation - z_CD4_exhaustion,
  NK_net_z   = z_NK_activation,
  MacDC_dysfunction_penalty_z = -z_frac_HAVCR2_of_MacDC,
  Function_z = rowMeans(cbind(CD8_net_z, CD4_net_z, NK_net_z, MacDC_dysfunction_penalty_z), na.rm = TRUE),
  Surveillance = Density_z * Function_z
)
write.csv(pm, file.path(out_dir, "STIC_immune_surveillance_score_joint_IFN_Metab.csv"), row.names = FALSE)

## 4. Primary contrast requested: double-High vs double-Low, paired by patient --
w <- pm %>% filter(joint_group %in% c("LowIFN_LowMetab","HighIFN_HighMetab")) %>%
  select(sample_id, joint_group, Density_z, Function_z, Surveillance) %>%
  pivot_wider(names_from = joint_group, values_from = c(Density_z, Function_z, Surveillance))
w <- w[complete.cases(w), ]
cat(sprintf("\n==== Double-High (IFN+Metab) vs Double-Low: paired Wilcoxon, n=%d patients ====\n", nrow(w)))
for (metric in c("Density_z", "Function_z", "Surveillance")) {
  lo <- w[[paste0(metric, "_LowIFN_LowMetab")]]; hi <- w[[paste0(metric, "_HighIFN_HighMetab")]]
  test <- wilcox.test(hi, lo, paired = TRUE)
  cat(sprintf("%-12s Double-Low mean=%.3f Double-High mean=%.3f | paired Wilcoxon p = %.4g\n",
              metric, mean(lo), mean(hi), test$p.value))
}

## 5. Faceted quadrant plot: one panel per joint group ------------------------
group_labels <- c(LowIFN_LowMetab = "IFN Low, Metab Low", LowIFN_HighMetab = "IFN Low, Metab High",
                   HighIFN_LowMetab = "IFN High, Metab Low", HighIFN_HighMetab = "IFN High, Metab High")
pm$joint_group_label <- factor(group_labels[as.character(pm$joint_group)], levels = group_labels)
pm$low_n_flag <- pm$n_cells < 5

p_facet <- ggplot(pm, aes(x = Density_z, y = Function_z, label = sample_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(color = Surveillance, shape = low_n_flag), size = 3) +
  ggrepel::geom_text_repel(size = 2.6, max.overlaps = 20) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
  scale_color_gradient2(low = "#2C7FB8", mid = "grey80", high = "#D7191C", midpoint = 0) +
  facet_wrap(~joint_group_label, nrow = 2) +
  theme_bw(base_size = 12) +
  labs(title = "Per-patient immune surveillance by joint IFN x Metabolic tertile state",
       subtitle = "X marks a patient-group cell with <5 cells (low-confidence estimate: STIC17, IFN High/Metab Low)",
       x = "Density_z (CD8+CD4+NK+MacDC infiltrate)",
       y = "Function_z (net activation - exhaustion - Mac/DC dysfunction)",
       color = "Surveillance\nscore")
ggsave(file.path(out_dir, "STIC_immune_surveillance_quadrant_joint_4panel.png"), p_facet, width = 11, height = 9.5, dpi = 300)

## 6. Bonus: single overlaid plot with each patient's 4-state trajectory ------
p_traj <- ggplot(pm, aes(x = Density_z, y = Function_z, group = sample_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_path(color = "grey50", alpha = 0.5, linewidth = 0.4,
            arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  geom_point(aes(color = joint_group_label, shape = low_n_flag), size = 2.6) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  labs(title = "Per-patient trajectory across joint IFN x Metabolic tertile states",
       subtitle = "Path follows Low_Low -> Low_High -> High_Low -> High_High for each patient",
       x = "Density_z", y = "Function_z", color = "Joint tertile state")
ggsave(file.path(out_dir, "STIC_immune_surveillance_quadrant_joint_trajectory.png"), p_traj, width = 8, height = 6.5, dpi = 300)

message("\nDone. Outputs written to: ", out_dir)
