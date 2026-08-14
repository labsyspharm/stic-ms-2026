# Step 3: merge all 16 samples and normalize with PFlog1pPF.
# Normalization must run AFTER merging -- the proportional-fitting target
# (mean total count) has to be computed across all pooled cells, otherwise


library(Seurat)
source("normalize_pflog1ppf.R")

samples <- readRDS("samples_list_filtered.rds")

seu <- merge(samples[[1]], samples[-1], add.cell.ids = names(samples))
seu <- JoinLayers(seu)

# Check whether counts are confounded with segmentation/cell area before
# trusting total-count-only depth normalization -- common in imaging-based
# spatial data. If |r| is large, consider area as an additional covariate
# (e.g. regress out during scaling, or use area instead of total counts as
# the PF target sum).
area_col <- intersect(c("cell_area", "nucleus_area"), colnames(seu@meta.data))
if (length(area_col) > 0) {
  r <- cor(seu$nCount_Xenium, seu[[area_col[1]]], use = "complete.obs")
  cat("Correlation(counts, ", area_col[1], ") = ", round(r, 3), "\n", sep = "")
}

counts <- GetAssayData(seu, assay = "Xenium", layer = "counts")
normalized <- pf_log1p_pf(counts)
seu <- SetAssayData(seu, assay = "Xenium", layer = "data", new.data = normalized)

saveRDS(seu, "merged_normalized.rds")

##### read the normalized sample file:
seu <- readRDS("merged_normalized.rds")

# pick a manageable slice: first 20 genes, first 100 cells
counts_slice <- GetAssayData(seu, assay = "Xenium", layer = "counts")[1:20, 1:100] #untouched raw transcripts counts
data_slice   <- GetAssayData(seu, assay = "Xenium", layer = "data")[1:20, 1:100] # PFlog1pPF- normalized values

# transpose so cells are rows, genes are columns -- more spreadsheet-like
View(as.matrix(t(counts_slice)))
View(as.matrix(t(data_slice)))






