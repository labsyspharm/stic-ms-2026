# PCA + Harmony + UMAP for the FULL ~6.39M-cell dataset, colored by the
# cluster labels already assigned via label transfer (predicted.id) --
# no FindNeighbors/FindClusters needed since those labels already exist.

library(Seurat)
library(harmony)

seu <- readRDS("full_dataset_with_transferred_clusters_lowclusters.rds")

seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = rownames(seu), npcs = 30)

seu <- RunHarmony(seu, group.by.vars = "run_id",
                   reduction.use = "pca", reduction.save = "harmony")

seu <- RunUMAP(seu, reduction = "harmony", dims = 1:10)

saveRDS(seu, "full_dataset_transferred_clusters_with_pca_umap.rds")
cat("Wrote full_dataset_transferred_clusters_with_pca_umap.rds\n")
