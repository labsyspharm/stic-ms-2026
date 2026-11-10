# Step 4 (PFlog1pPF variant): same pipeline as 04_cluster.R, but clusters on
# the "data" layer (PFlog1pPF, from 03_merge_normalize.R) instead of
# "data_od" (overdispersion/CLR). merged_normalized_overdispersion.rds
# carries both layers, so no materialization step is needed here -- "data"
# is already the finished PFlog1pPF result.

library(Seurat)
library(harmony)

seu <- readRDS("merged_normalized_overdispersion.rds")

# 5000 cells PER SAMPLE this time (16 samples x 5000 = up to 80,000 cells
# total), rather than 5000 pooled across the whole dataset -- guarantees
# every sample is represented even though sample sizes vary a lot
set.seed(1)
seu <- seu[, unlist(lapply(split(colnames(seu), seu$sample_id), function(cells) {
  sample(cells, min(5000, length(cells)))
}))]

seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = rownames(seu), npcs = 30)

seu <- RunHarmony(seu, group.by.vars = "run_id",
                   reduction.use = "pca", reduction.save = "harmony")

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:10)
seu <- FindClusters(seu, resolution = 0.8)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:10)

saveRDS(seu, "clustered_subsample5000persample_pflog.rds")

seu <- readRDS("clustered_subsample5000_pflog.rds")

data <- GetAssayData(seu, assay = "Xenium", layer = "data")


seu <- readRDS("merged_normalized.rds")  # adjust path/filename as needed

data <- GetAssayData(seu, assay = "Xenium", layer = "data")

# pull a small slice and convert to a plain matrix -- View() gives you a
# proper spreadsheet-style window instead of console text
mat <- as.matrix(data[1:20, 1:50])   # first 20 genes, first 50 cells
View(mat)


# --- Lower-resolution re-clustering ------------------------------------
# FindNeighbors already computed the harmony-based neighbor graph inside
# clustered_subsample5000persample_pflog.rds -- only FindClusters needs to
# rerun for a different resolution, no need to redo ScaleData/PCA/Harmony/
# FindNeighbors/RunUMAP (UMAP coordinates don't depend on the resolution
# used for cluster labels, just the harmony embedding, which is unchanged).
seu_lowres <- readRDS("clustered_subsample5000persample_pflog.rds")
seu_lowres <- FindClusters(seu_lowres, resolution = 0.4)
cat("Number of clusters at resolution 0.4:",
    length(unique(seu_lowres$seurat_clusters)), "\n")
saveRDS(seu_lowres, "clustered_subsample5000persample_pflog_lowclusters.rds")






 
