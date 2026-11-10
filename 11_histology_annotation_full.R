# Integrate histopathological lesion annotations (coordinate[s]_annotation_
# STIC.csv, Xenium Explorer polygon exports) with the full dataset.
#
# Each of the 12 raw runs has one of these files, giving one or more named
# lesion polygons (e.g. "STIC-BLAD-1", "STIC-FLAT", "STIL-2", "p53
# signatures-3", "Tumor"). Per-focus labels are kept distinct (per user's
# choice -- "STIC-BLAD-1" and "STIC-BLAD-2" are NOT collapsed to one label).
#
# Rule: a cell's centroid (x_centroid, y_centroid, already carried in
# meta.data since 01_load_and_split.R) inside exactly one polygon gets that
# polygon's label; inside zero polygons gets "Normal"; inside more than one
# (overlap) is flagged ambiguous and left NA.
#
# Point-in-polygon is done per run_id (not per final sample_id), since the
# annotation coordinates are in the same raw per-slide coordinate space
# that a run's two specimens (if split) both live in.
#
# Plots use scale_y_reverse() because Xenium centroid coordinates are in
# image-space (origin top-left, y increasing downward, same as Xenium
# Explorer / the physical slide) -- ggplot's default y-increasing-upward
# would otherwise render the tissue vertically flipped.

library(Seurat)
library(sp)
library(ggplot2)
library(readr)
library(dplyr)

out_dir <- "histology_annotation"
plot_dir <- file.path(out_dir, "plots")

find_annotation_file <- function(run_id) {
  candidates <- c(
    file.path("rawdata", run_id, "coordinate_annotation_STIC.csv"),
    file.path("rawdata", run_id, "coordinates_annotation_STIC.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) stop("No annotation file found for run ", run_id)
  hit[1]
}

assign_histology <- function(x_centroid, y_centroid, poly) {
  sel_names <- unique(poly$Selection)
  inside <- vapply(sel_names, function(s) {
    p <- poly[poly$Selection == s, ]
    sp::point.in.polygon(x_centroid, y_centroid, p$X, p$Y) > 0
  }, logical(length(x_centroid)))
  if (is.null(dim(inside))) inside <- matrix(inside, ncol = length(sel_names))
  n_hits <- rowSums(inside)

  label <- rep("Normal", length(x_centroid))
  single_hit <- n_hits == 1
  label[single_hit] <- sel_names[max.col(inside[single_hit, , drop = FALSE])]
  label[n_hits > 1] <- NA_character_

  list(label = label, n_ambiguous = sum(n_hits > 1))
}

cat("Loading full dataset...\n")
seu <- readRDS("full_dataset_transferred_clusters_with_pca_umap.rds")

meta <- seu@meta.data
run_ids <- sort(unique(meta$run_id))
cat("Found", length(run_ids), "runs:", paste(run_ids, collapse = ", "), "\n\n")

histology_col <- setNames(rep(NA_character_, nrow(meta)), rownames(meta))
summary_rows <- list()

for (run_id in run_ids) {
  idx <- which(meta$run_id == run_id)
  cell_ids <- rownames(meta)[idx]

  ann_path <- find_annotation_file(run_id)
  poly <- read.csv(ann_path, comment.char = "#")

  res <- assign_histology(meta$x_centroid[idx], meta$y_centroid[idx], poly)
  histology_col[cell_ids] <- res$label

  tab <- table(res$label, useNA = "ifany")
  cat("===", run_id, "(", ann_path, ") ===\n")
  cat("n cells:", length(idx), " | ambiguous:", res$n_ambiguous, "\n")
  print(tab)
  cat("\n")

  summary_rows[[run_id]] <- data.frame(
    run_id = run_id, label = names(tab), n_cells = as.integer(tab)
  )

  # per-run plot, in tissue orientation
  df <- data.frame(x_centroid = meta$x_centroid[idx],
                    y_centroid = meta$y_centroid[idx],
                    histology = res$label)
  p <- ggplot(df, aes(x = x_centroid, y = y_centroid, color = histology)) +
    geom_point(size = 0.3) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = paste("Histology annotation --", run_id)) +
    theme_minimal()

  png(file.path(plot_dir, paste0(run_id, "_histology.png")),
      width = 1400, height = 1100, res = 150)
  print(p)
  dev.off()
}

seu$histology <- histology_col[colnames(seu)]

cat("\nOverall label counts (all runs combined):\n")
print(table(seu$histology, useNA = "ifany"))

saveRDS(seu, "full_dataset_with_histology_annotation.rds")
cat("\nWrote full_dataset_with_histology_annotation.rds\n")

summary_df <- bind_rows(summary_rows)
write_csv(summary_df, file.path(out_dir, "histology_label_counts_per_run.csv"))
cat("Wrote", file.path(out_dir, "histology_label_counts_per_run.csv"), "\n")
cat("Wrote per-run plots to", plot_dir, "\n")
