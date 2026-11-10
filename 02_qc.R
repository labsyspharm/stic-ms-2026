# Step 2: per-sample QC. Run BEFORE merging -- depth, segmentation quality,
# and negative-control rates can differ by run, so thresholds are assessed
# per sample first.

library(Seurat)
library(dplyr)

samples <- readRDS("samples_list_raw.rds")

qc_table <- bind_rows(lapply(names(samples), function(id) {
  s <- samples[[id]]
  data.frame(
    sample_id       = id,
    run_id          = unique(s$run_id),
    n_cells         = ncol(s),
    median_counts   = median(s$nCount_Xenium),
    median_genes    = median(s$nFeature_Xenium),
    median_ctrl     = median(s$nCount_controls),
    # fraction of a cell's total signal coming from negative controls --
    # should be near zero; a run with elevated values has a probe/assay issue
    median_ctrl_frac = median(s$nCount_controls /
                                 (s$nCount_Xenium + s$nCount_controls),
                               na.rm = TRUE)
  )
}))
print(qc_table)
write.csv(qc_table, "qc_summary_by_sample.csv", row.names = FALSE)

# Adjust these after reviewing qc_table -- Xenium panels run lower depth
# than scRNA-seq, so thresholds are gentler than typical single-cell QC.
MIN_COUNTS <- 10
MIN_FEATURES <- 5
MAX_CTRL_FRAC <- 0.1

samples_filtered <- lapply(samples, function(s) {
  ctrl_frac <- s$nCount_controls / (s$nCount_Xenium + s$nCount_controls)
  keep <- s$nCount_Xenium >= MIN_COUNTS &
    s$nFeature_Xenium >= MIN_FEATURES &
    ctrl_frac <= MAX_CTRL_FRAC &
    # nucleus_count > 1 means the segmentation boundary caught more than one
    # nucleus -- these cells have several-fold more counts/genes/area than
    # single-nucleus cells, consistent with multiple real cells merged into
    # one polygon rather than a genuine multinucleate cell
    s$nucleus_count <= 1
  s[, keep]
})

cat("Cells retained per sample:\n")
print(sapply(samples_filtered, ncol) / sapply(samples, ncol))

saveRDS(samples_filtered, "samples_list_filtered.rds")

samples_filtered <- readRDS("samples_list_filtered.rds")
