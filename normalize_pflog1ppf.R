# PFlog1pPF normalization (Booeshaghi, Hallgrimsdottir & Pachter 2022,
# "Normalization for sampled count data", biorxiv 2022.05.06.490859)
#
# Proportional Fitting -> log1p -> Proportional Fitting
# Removes depth-correlation left behind by plain log1p(CPM), without the
# false-discovery inflation seen with sctransform/Pearson residuals.

proportional_fitting <- function(mat, target_sum = NULL) {
  # mat: genes x cells (dgCMatrix or matrix)
  cell_sums <- Matrix::colSums(mat)
  if (is.null(target_sum)) target_sum <- mean(cell_sums)
  sf <- target_sum / cell_sums
  mat_pf <- mat %*% Matrix::Diagonal(x = sf)
  dimnames(mat_pf) <- dimnames(mat)
  mat_pf
}

pf_log1p_pf <- function(mat) {
  mat |> proportional_fitting() |> log1p() |> proportional_fitting()
}

# --- Usage with a Seurat object ---
# counts <- GetAssayData(seu, assay = "Xenium", layer = "counts")
# normalized <- pf_log1p_pf(counts)
# seu <- SetAssayData(seu, assay = "Xenium", layer = "data", new.data = normalized)
