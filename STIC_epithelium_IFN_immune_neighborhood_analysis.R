## STIC epithelium: overall IFN score (chronic+acute combined) vs spatial immune neighborhood
##
## Index cells:   epithelial (predicted.id %in% c("4","7")), histology starts "STIC-",
##                excluding samples STIC12/STIC16.
## Neighbor pool: ALL cells (any type/histology) within the same 14 retained samples
##                (a true spatial neighborhood, not restricted to a cell type).
## KNN:           k=20 nearest neighbors by (x_centroid, y_centroid) Euclidean distance,
##                computed PER SAMPLE (coordinates are only meaningful within one section)
##                via RANN::nn2 (compiled KD-tree). Self is excluded after a distance~0
##                sanity check (the query cell is present in its own sample's candidate
##                pool, so it should always be its own nearest neighbor at distance 0).
## Immune gates:  computed on normalized/log-transformed expression ("data" layer), NOT
##                raw counts -- a raw-count magnitude threshold would be biased by
##                per-cell library size/cell size; log-normalized data corrects for that.
##                Presence/absence (>0) itself is unaffected by normalization since it's
##                a positive monotonic rescale, but we standardize on the normalized layer
##                throughout for consistency with the UCell functional scores below.
##
## Signatures (from STIC_priority_gene_panel_v5.xlsx):
##   Overall IFN (chronic+acute combined): STAT1,IRF1,IRF2,SOCS1,SOCS3,ISG15,MX1,OAS1,IFIT2,IFIT3,IFITM3,LY6E
##   CD8+ T cell lineage:   CD8A or CD8B (marker-only, no CD3 requirement)
##   CD4+ T cell lineage:   CD4 (marker-only, no CD3 requirement)
##   CD8+ exhaustion:       UCell(CTLA4, HAVCR2, LAG3, TIGIT) within CD8+ neighbors
##   CD8+ activation:       UCell(GNLY, NKG7, GZMB, PRF1) within CD8+ neighbors
##   CD4+ exhaustion:       UCell(CTLA4, HAVCR2, LAG3, TIGIT) within CD4+ neighbors
##   CD4+ activation:       UCell(GNLY, NKG7, GZMB, PRF1) within CD4+ neighbors
##   Macrophage/DC lineage: CD68 or CD163 or CD14 or ITGAM or ITGAX
##   HAVCR2+ Mac/DC:        Macrophage/DC lineage AND HAVCR2+
##   NK identity:           CD3-, CD8-, CD4- AND (KLRC1 or KLRK1)
##   NK activation:         UCell(GNLY, NKG7, GZMB, PRF1) within NK+ neighbors

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(RANN)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

data_dir <- "c:/Users/tak415/OneDrive - Harvard University/Documents/CyCIF_Geomx_analysis_LSP_HMS/CYCIF_ORION_atlas_2.0/Xenium_analysis"
out_dir  <- data_dir
set.seed(1)

excluded_samples <- c("STIC12", "STIC16")

chronic_genes      <- c("STAT1", "IRF1", "IRF2", "SOCS1", "SOCS3")
acute_genes        <- c("ISG15", "MX1", "OAS1", "IFIT2", "IFIT3", "IFITM3", "LY6E")
overall_ifn_genes  <- c(chronic_genes, acute_genes)
exhaustion_genes   <- c("CTLA4", "HAVCR2", "LAG3", "TIGIT")
activation_genes   <- c("GNLY", "NKG7", "GZMB", "PRF1")
k_neighbors <- 20

## 1. Load -------------------------------------------------------------------
message("Loading Seurat object...")
obj <- readRDS(file.path(data_dir, "full_dataset_with_histology_annotation.rds"))
DefaultAssay(obj) <- "Xenium"
obj <- subset(obj, subset = !(sample_id %in% excluded_samples))
obj <- NormalizeData(obj, verbose = FALSE)

need_genes <- unique(c(overall_ifn_genes, exhaustion_genes, activation_genes,
                        "CD3D","CD3E","CD4","CD8A","CD8B","CD68","CD163","CD14","ITGAM","ITGAX",
                        "HAVCR2","KLRC1","KLRK1"))
stopifnot(all(need_genes %in% rownames(obj)))

## 2. Immune lineage gates (normalized data > 0) + functional UCell scores ---
## computed once for ALL retained cells (candidate neighbor pool), any cell type/histology.
message("Computing immune lineage gates and functional scores for full candidate pool...")
expr <- GetAssayData(obj, assay = "Xenium", layer = "data")
pos <- function(gene) expr[gene, ] > 0

md <- obj@meta.data
md$sample_id <- as.character(md$sample_id)
md$CD8_pos   <- pos("CD8A") | pos("CD8B")
md$CD4_pos   <- pos("CD4")
md$NK_pos    <- !(pos("CD3D") | pos("CD3E")) & !md$CD8_pos & !md$CD4_pos & (pos("KLRC1") | pos("KLRK1"))
md$MacDC_pos <- pos("CD68") | pos("CD163") | pos("CD14") | pos("ITGAM") | pos("ITGAX")
md$HAVCR2_pos <- pos("HAVCR2")

cat("\nCandidate pool cells:", nrow(md), "across", length(unique(md$sample_id)), "samples\n")
cat("CD8+ cells:", sum(md$CD8_pos), "| CD4+ cells:", sum(md$CD4_pos), "| NK+ cells:", sum(md$NK_pos),
    "| Mac/DC+ cells:", sum(md$MacDC_pos), "| HAVCR2+ of Mac/DC:", sum(md$MacDC_pos & md$HAVCR2_pos), "\n")

## Exhaustion/ActivationFn scores are only ever read out for CD8+/CD4+/NK+ neighbor
## cells downstream (everything else is masked to NA at aggregation time), so
## restrict the UCell computation to that subset instead of all ~5M candidate
## cells -- same output, much less compute.
immune_mask <- md$CD8_pos | md$CD4_pos | md$NK_pos
cat("\nCD8+ or CD4+ or NK+ cells (subset for functional UCell scoring):", sum(immune_mask), "/", nrow(md), "\n")

obj_immune <- subset(obj, cells = rownames(md)[immune_mask])
obj_immune <- AddModuleScore_UCell(
  obj_immune, features = list(Exhaustion = exhaustion_genes, ActivationFn = activation_genes),
  assay = "Xenium", slot = "data", name = ""
)
## Match back by cell name explicitly (not by row position) -- subset() is not
## guaranteed to preserve row order, and a silent misalignment here would
## silently corrupt every downstream exhaustion/activation readout.
match_idx <- match(rownames(md), rownames(obj_immune@meta.data))
md$Exhaustion   <- obj_immune@meta.data$Exhaustion[match_idx]
md$ActivationFn <- obj_immune@meta.data$ActivationFn[match_idx]
stopifnot(identical(is.na(md$Exhaustion), !immune_mask))  # NA exactly where not CD8+/CD4+/NK+
rm(obj_immune); gc(verbose = FALSE)

## 3. Overall IFN score for STIC epithelial index cells -----------------------
is_epithelial <- md$predicted.id %in% c("4", "7")
is_stic <- grepl("^STIC-", md$histology) & !is.na(md$histology)
idx_mask <- is_epithelial & is_stic

message(sprintf("STIC epithelial index cells: %d", sum(idx_mask)))
obj_idx <- subset(obj, cells = rownames(md)[idx_mask])
obj_idx <- AddModuleScore_UCell(
  obj_idx, features = list(OverallIFN = overall_ifn_genes),
  assay = "Xenium", slot = "data", name = ""
)

idx_scores <- data.frame(
  cell = rownames(obj_idx@meta.data),
  sample_id = as.character(obj_idx@meta.data$sample_id),
  x = obj_idx@meta.data$x_centroid,
  y = obj_idx@meta.data$y_centroid,
  OverallIFN = obj_idx@meta.data$OverallIFN,
  stringsAsFactors = FALSE
)

q_tert <- quantile(idx_scores$OverallIFN, probs = c(1/3, 2/3))
idx_scores$IFN_tertile <- cut(idx_scores$OverallIFN, breaks = c(-Inf, q_tert[1], q_tert[2], Inf),
                               labels = c("Low", "Mid", "High"))

cat("\nSTIC epithelial index cells:", nrow(idx_scores), "\n")
cat("Tertile counts:\n"); print(table(idx_scores$IFN_tertile))
cat("Tertile cutoffs (OverallIFN UCell score):\n"); print(q_tert)
cat("Samples represented among index cells:\n"); print(table(idx_scores$sample_id))

rm(obj_idx); gc(verbose = FALSE)

## 4. Per-sample KNN (k=20), candidate pool = ALL retained cells --------------
message("Running per-sample spatial KNN (RANN::nn2, KD-tree)...")
samples <- sort(unique(idx_scores$sample_id))

CD8_pos_vec      <- md$CD8_pos
CD4_pos_vec      <- md$CD4_pos
NK_pos_vec       <- md$NK_pos
MacDC_pos_vec    <- md$MacDC_pos
HAVCR2_MacDC_vec <- md$MacDC_pos & md$HAVCR2_pos
Exhaustion_vec   <- md$Exhaustion
ActivationFn_vec <- md$ActivationFn

result_list <- vector("list", length(samples))
names(result_list) <- samples

for (s in samples) {
  cand_mask   <- md$sample_id == s
  cand_coords <- as.matrix(md[cand_mask, c("x_centroid", "y_centroid")])
  cand_rows   <- which(cand_mask)  # global row indices into md

  q_mask   <- idx_scores$sample_id == s
  q_coords <- as.matrix(idx_scores[q_mask, c("x", "y")])
  q_cells  <- idx_scores$cell[q_mask]

  nn <- RANN::nn2(data = cand_coords, query = q_coords, k = k_neighbors + 1)

  self_ok <- nn$nn.dists[, 1] < 1e-6
  if (mean(self_ok) < 0.99) {
    warning(sprintf("Sample %s: only %.1f%% of query cells matched self as nearest neighbor (dist~0) -- check for duplicate coordinates.",
                     s, 100 * mean(self_ok)))
  }
  neighbor_idx_local  <- nn$nn.idx[, 2:(k_neighbors + 1), drop = FALSE]
  neighbor_idx_global <- matrix(cand_rows[neighbor_idx_local], nrow = nrow(neighbor_idx_local))

  CD8_mat   <- matrix(CD8_pos_vec[neighbor_idx_global],      nrow = nrow(neighbor_idx_global))
  CD4_mat   <- matrix(CD4_pos_vec[neighbor_idx_global],      nrow = nrow(neighbor_idx_global))
  NK_mat    <- matrix(NK_pos_vec[neighbor_idx_global],       nrow = nrow(neighbor_idx_global))
  MacDC_mat <- matrix(MacDC_pos_vec[neighbor_idx_global],    nrow = nrow(neighbor_idx_global))
  HAVCR2_mat<- matrix(HAVCR2_MacDC_vec[neighbor_idx_global], nrow = nrow(neighbor_idx_global))
  Exh_mat   <- matrix(Exhaustion_vec[neighbor_idx_global],   nrow = nrow(neighbor_idx_global))
  Act_mat   <- matrix(ActivationFn_vec[neighbor_idx_global], nrow = nrow(neighbor_idx_global))

  n_cd8   <- rowSums(CD8_mat)
  n_cd4   <- rowSums(CD4_mat)
  n_nk    <- rowSums(NK_mat)
  n_macdc <- rowSums(MacDC_mat)
  n_havcr2<- rowSums(HAVCR2_mat)

  Exh_cd8_mat <- Exh_mat; Exh_cd8_mat[!CD8_mat] <- NA
  Act_cd8_mat <- Act_mat; Act_cd8_mat[!CD8_mat] <- NA
  Exh_cd4_mat <- Exh_mat; Exh_cd4_mat[!CD4_mat] <- NA
  Act_cd4_mat <- Act_mat; Act_cd4_mat[!CD4_mat] <- NA
  Act_nk_mat  <- Act_mat; Act_nk_mat[!NK_mat]   <- NA

  cd8_exhaustion <- rowMeans(Exh_cd8_mat, na.rm = TRUE); cd8_exhaustion[n_cd8 == 0] <- NA
  cd8_activation <- rowMeans(Act_cd8_mat, na.rm = TRUE); cd8_activation[n_cd8 == 0] <- NA
  cd4_exhaustion <- rowMeans(Exh_cd4_mat, na.rm = TRUE); cd4_exhaustion[n_cd4 == 0] <- NA
  cd4_activation <- rowMeans(Act_cd4_mat, na.rm = TRUE); cd4_activation[n_cd4 == 0] <- NA
  nk_activation  <- rowMeans(Act_nk_mat,  na.rm = TRUE); nk_activation[n_nk == 0]   <- NA

  result_list[[s]] <- data.frame(
    cell = q_cells,
    sample_id = s,
    frac_CD8   = n_cd8   / k_neighbors,
    frac_CD4   = n_cd4   / k_neighbors,
    frac_NK    = n_nk    / k_neighbors,
    frac_MacDC = n_macdc / k_neighbors,
    frac_HAVCR2_of_all   = n_havcr2 / k_neighbors,
    frac_HAVCR2_of_MacDC = ifelse(n_macdc > 0, n_havcr2 / n_macdc, NA_real_),
    CD8_exhaustion = cd8_exhaustion,
    CD8_activation = cd8_activation,
    CD4_exhaustion = cd4_exhaustion,
    CD4_activation = cd4_activation,
    NK_activation  = nk_activation,
    stringsAsFactors = FALSE
  )
}

readouts <- bind_rows(result_list)
readouts <- idx_scores %>%
  select(cell, sample_id, OverallIFN, IFN_tertile) %>%
  inner_join(readouts, by = c("cell", "sample_id"))

cat("\nReadouts table dim:\n"); print(dim(readouts))
write.csv(readouts, file.path(out_dir, "STIC_epithelial_neighborhood_readouts.csv"), row.names = FALSE)
saveRDS(readouts, file.path(out_dir, "STIC_epithelial_neighborhood_readouts.rds"))

## 5. Stats -------------------------------------------------------------------
readout_cols <- c("frac_CD8", "frac_CD4", "frac_NK", "frac_MacDC", "frac_HAVCR2_of_all",
                   "frac_HAVCR2_of_MacDC", "CD8_exhaustion", "CD8_activation",
                   "CD4_exhaustion", "CD4_activation", "NK_activation")

fit_continuous <- function(col) {
  d <- readouts[!is.na(readouts[[col]]), c(col, "OverallIFN", "sample_id")]
  f <- as.formula(paste0(col, " ~ OverallIFN + (1|sample_id)"))
  m <- lmer(f, data = d)
  s <- summary(m)
  ci <- tryCatch(confint(m, parm = "OverallIFN", method = "Wald"), error = function(e) matrix(NA, 1, 2))
  coef_est <- s$coefficients["OverallIFN", "Estimate"]
  raw_cor <- cor(d$OverallIFN, d[[col]], method = "spearman")
  same_sign <- sign(raw_cor) == sign(coef_est)
  sign_ok <- same_sign || abs(raw_cor) <= 1e-6
  cat(sprintf("[sign check] %s: raw Spearman rho = %.4f | LMM coef = %.5f | signs agree: %s\n",
              col, raw_cor, coef_est, sign_ok))
  list(model = m, summary = s, ci = ci, n = nrow(d), raw_cor = raw_cor, sign_ok = sign_ok)
}

fit_highlow <- function(col) {
  d <- readouts %>% filter(IFN_tertile %in% c("Low", "High"), !is.na(.data[[col]]))
  d$IFN_tertile <- factor(d$IFN_tertile, levels = c("Low", "High"))
  stopifnot(identical(levels(d$IFN_tertile), c("Low", "High")))

  f <- as.formula(paste0(col, " ~ IFN_tertile + (1|sample_id)"))
  m <- lmer(f, data = d)
  s <- summary(m)
  raw_means <- tapply(d[[col]], d$IFN_tertile, mean)
  raw_diff_pooled <- unname(raw_means["High"] - raw_means["Low"])
  coef_est <- s$coefficients["IFN_tertileHigh", "Estimate"]

  sm <- d %>% group_by(sample_id, IFN_tertile) %>%
    summarise(mean_val = mean(.data[[col]]), n = n(), .groups = "drop")
  sm_wide <- sm %>% select(sample_id, IFN_tertile, mean_val) %>%
    tidyr::pivot_wider(names_from = IFN_tertile, values_from = mean_val)
  sm_wide <- sm_wide[complete.cases(sm_wide), ]
  paired_test <- if (nrow(sm_wide) >= 2) wilcox.test(sm_wide$High, sm_wide$Low, paired = TRUE) else NULL
  patient_weighted_diff <- if (nrow(sm_wide) >= 1) mean(sm_wide$High - sm_wide$Low) else NA_real_

  sign_ok_pooled <- (sign(raw_diff_pooled) == sign(coef_est)) || abs(raw_diff_pooled) <= 1e-8
  sign_ok_patient <- (sign(patient_weighted_diff) == sign(coef_est)) || abs(patient_weighted_diff) <= 1e-8
  cat(sprintf("[sign check] %s (High vs Low): pooled raw diff = %.5f | patient-weighted diff = %.5f | LMM coef = %.5f | pooled agree: %s | patient-weighted agree: %s\n",
              col, raw_diff_pooled, patient_weighted_diff, coef_est, sign_ok_pooled, sign_ok_patient))

  list(model = m, summary = s, raw_means = raw_means, sm_wide = sm_wide, paired_test = paired_test,
       raw_diff_pooled = raw_diff_pooled, patient_weighted_diff = patient_weighted_diff,
       sign_ok_pooled = sign_ok_pooled, sign_ok_patient = sign_ok_patient)
}

cat("\n\n================ CONTINUOUS MODELS (readout ~ OverallIFN + (1|sample_id)) ================\n")
continuous_fits <- list()
for (col in readout_cols) {
  cat("\n----", col, "----\n")
  fit <- tryCatch(fit_continuous(col), error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(fit)) {
    continuous_fits[[col]] <- fit
    print(fit$summary$coefficients)
    cat("95% CI:\n"); print(fit$ci)
    cat("n =", fit$n, "\n")
  }
}

cat("\n\n================ HIGH vs LOW TERTILE MODELS ================\n")
highlow_fits <- list()
for (col in readout_cols) {
  cat("\n----", col, "----\n")
  fit <- tryCatch(fit_highlow(col), error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(fit)) {
    highlow_fits[[col]] <- fit
    cat("Raw means:\n"); print(fit$raw_means)
    print(fit$summary$coefficients)
    if (!is.null(fit$paired_test)) {
      cat(sprintf("Per-patient PAIRED Wilcoxon (n=%d patients): p = %.4g\n",
                   nrow(fit$sm_wide), fit$paired_test$p.value))
    } else {
      cat("Per-patient paired Wilcoxon: not enough patients with both tertiles present.\n")
    }
  }
}

## 6. Stats summary CSV --------------------------------------------------------
summary_rows <- lapply(readout_cols, function(col) {
  cf <- continuous_fits[[col]]; hf <- highlow_fits[[col]]
  genuinely_unresolved <- !is.null(hf) && !hf$sign_ok_pooled && !hf$sign_ok_patient
  pooled_vs_patient_disagree <- !is.null(hf) && !hf$sign_ok_pooled && hf$sign_ok_patient
  patient_confirmed <- !is.null(hf) && !is.null(cf) && !is.null(hf$paired_test) &&
    hf$paired_test$p.value < 0.05 && cf$summary$coefficients["OverallIFN","Pr(>|t|)"] < 0.05 &&
    hf$summary$coefficients["IFN_tertileHigh","Pr(>|t|)"] < 0.05
  interp <- if (is.null(cf) || is.null(hf)) {
    "ERROR: model fit failed"
  } else if (genuinely_unresolved) {
    "UNRESOLVED: direction disagrees even after equal patient weighting"
  } else if (pooled_vs_patient_disagree && patient_confirmed) {
    "Patient-level confirmed despite pooled cell-level sign disagreement (pseudo-replication imbalance)"
  } else if (patient_confirmed) {
    "Robust: significant at cell level (both LMMs) AND patient level (paired Wilcoxon)"
  } else if (!is.null(cf) && !is.null(hf) &&
             cf$summary$coefficients["OverallIFN","Pr(>|t|)"] < 0.05 &&
             hf$summary$coefficients["IFN_tertileHigh","Pr(>|t|)"] < 0.05) {
    "Suggestive: significant at cell level, not confirmed at patient level (n=14)"
  } else {
    "No significant association"
  }
  data.frame(
    readout = col,
    continuous_lmm_p = if (!is.null(cf)) cf$summary$coefficients["OverallIFN","Pr(>|t|)"] else NA,
    highlow_lmm_p = if (!is.null(hf)) hf$summary$coefficients["IFN_tertileHigh","Pr(>|t|)"] else NA,
    paired_wilcoxon_p = if (!is.null(hf) && !is.null(hf$paired_test)) hf$paired_test$p.value else NA,
    n_patients_paired = if (!is.null(hf)) nrow(hf$sm_wide) else NA,
    interpretation = interp,
    stringsAsFactors = FALSE
  )
})
summary_df <- bind_rows(summary_rows)
write.csv(summary_df, file.path(out_dir, "STIC_epithelial_neighborhood_stats_summary.csv"), row.names = FALSE)
cat("\n\nStats summary:\n"); print(summary_df)

## 7. Two-panel plots (cell-level + patient-level) directly ------------------
readout_labels <- c(
  frac_CD8 = "Fraction of 20 neighbors: CD8+",
  frac_CD4 = "Fraction of 20 neighbors: CD4+",
  frac_NK = "Fraction of 20 neighbors: NK+",
  frac_MacDC = "Fraction of 20 neighbors: Mac/DC+",
  frac_HAVCR2_of_all = "Fraction of 20 neighbors: HAVCR2+ Mac/DC",
  frac_HAVCR2_of_MacDC = "Fraction of Mac/DC neighbors: HAVCR2+",
  CD8_exhaustion = "Mean CD8+ exhaustion UCell score (neighbors)",
  CD8_activation = "Mean CD8+ activation UCell score (neighbors)",
  CD4_exhaustion = "Mean CD4+ exhaustion UCell score (neighbors)",
  CD4_activation = "Mean CD4+ activation UCell score (neighbors)",
  NK_activation = "Mean NK activation UCell score (neighbors)"
)

make_twopanel <- function(col, ylab) {
  d <- readouts %>% filter(IFN_tertile %in% c("Low", "High"), !is.na(.data[[col]]))
  d$IFN_tertile <- factor(d$IFN_tertile, levels = c("Low", "High"))

  pm <- d %>% group_by(sample_id, IFN_tertile) %>%
    summarise(mean_val = mean(.data[[col]]), .groups = "drop")
  pm_wide <- pm %>% tidyr::pivot_wider(names_from = IFN_tertile, values_from = mean_val)
  pm_wide <- pm_wide[complete.cases(pm_wide), ]
  paired_test <- wilcox.test(pm_wide$High, pm_wide$Low, paired = TRUE)

  dir_lookup <- pm_wide %>% mutate(direction = ifelse(High < Low, "High < Low", "High >= Low"))
  pm <- pm %>% left_join(dir_lookup %>% select(sample_id, direction), by = "sample_id")

  p_left <- ggplot(d, aes(x = IFN_tertile, y = .data[[col]], color = IFN_tertile)) +
    geom_jitter(width = 0.15, alpha = 0.12, size = 0.6, shape = 16) +
    geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", color = "black", alpha = 0.6) +
    scale_color_manual(values = c(Low = "#F8766D", High = "#00BFC4")) +
    theme_bw(base_size = 12) +
    labs(title = "Cell level", x = "Index cell's own IFN tertile", y = ylab) +
    theme(legend.position = "none")

  p_right <- ggplot(pm, aes(x = IFN_tertile, y = mean_val, group = sample_id, color = direction)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c("High < Low" = "grey30", "High >= Low" = "#D55E00")) +
    theme_bw(base_size = 12) +
    labs(title = sprintf("Patient level (n=%d)", nrow(pm_wide)),
         x = "IFN tertile group averaged within patient",
         y = paste0("Per-patient mean ", col, " (neighbors)"),
         color = NULL) +
    theme(legend.position = "bottom")

  combined <- p_left + p_right +
    plot_annotation(
      title = paste0(col, ": High vs Low overall IFN (STIC epithelium neighbors)"),
      subtitle = sprintf("Per-patient paired Wilcoxon p = %s", signif(paired_test$p.value, 3))
    )

  ggsave(file.path(out_dir, paste0("STIC_neighborhood_", col, "_highlow_twopanel.png")),
         combined, width = 9.5, height = 5, dpi = 300)
}

for (col in readout_cols) make_twopanel(col, readout_labels[[col]])

message("\nDone. Outputs written to: ", out_dir)
