# Step 5 (label transfer): project the 5000-cell subsample's cluster labels
# onto the full ~6.39M-cell dataset, rather than redoing Louvain clustering
# from scratch at full scale. Uses Seurat's standard reference-mapping
# workflow: project every cell into the subsample's PCA space (same gene
# loadings/scaling the subsample was clustered in), find anchors against
# the labeled reference cells, and assign each cell the best-matching
# cluster plus a confidence score.

library(Seurat)

# default future.globals.maxSize (500MB) is too small for a ~6.4M-cell
# query; this machine has ~128GB RAM so raise it generously
options(future.globals.maxSize = 20 * 1024^3)

ref <- readRDS("clustered_subsample5000persample_pflog_lowclusters.rds")

# full dataset, same normalization (PFlog1pPF) the reference was clustered
# on -- deliberately NOT excluding the 80,000 reference cells from this;
# the small redundant compute cost buys a free sanity check below
query <- readRDS("normalisation/merged_normalized.rds")

anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  dims = 1:10,
  reference.reduction = "pca",
  features = rownames(ref)
)

predictions <- TransferData(
  anchorset = anchors,
  refdata = ref$seurat_clusters,
  dims = 1:10
)

query <- AddMetaData(query, predictions)

cat("Predicted cluster sizes across full dataset:\n")
print(table(query$predicted.id))

cat("\nPrediction confidence (prediction.score.max) summary:\n")
print(summary(query$prediction.score.max))

# sanity check: the reference cells are present in query too (not
# excluded) -- their predicted label should match their known original
# cluster almost perfectly if the projection is working correctly
ref_cells_in_query <- query@meta.data[colnames(ref), ]
agreement <- mean(as.character(ref_cells_in_query$predicted.id) ==
                     as.character(ref$seurat_clusters))
cat("\nSanity check -- agreement between predicted.id and known\n")
cat("seurat_clusters, restricted to the 80,000 reference cells:", agreement, "\n")

saveRDS(query, "full_dataset_with_transferred_clusters_lowclusters.rds")
cat("\nWrote full_dataset_with_transferred_clusters_lowclusters.rds\n")
