# Step 1: load each of the 12 Xenium runs and resolve them down to 16
# per-sample Seurat objects. Fill in manifest_template.csv (rename it
# manifest.csv) with one row per run before running this.
#
# split_method per row:
#   "single"        - run has exactly one sample, no splitting needed
#   "annotation"     - Xenium Explorer "Export Selections" (coordinates.csv):
#                      polygon vertices per named selection. Cells are
#                      assigned to whichever polygon contains their centroid.
#   "spatial_kmeans" - no annotation exists; split 2 tissue pieces by x/y
#                      centroid clustering (inspect the plot it produces!)

library(Seurat)
library(dplyr)
library(sp)

CONTROL_PREFIXES <- c("NegControlProbe_", "NegControlCodeword_",
                       "BLANK_", "UnassignedCodeword_", "antisense_")

# Xenium's feature-barcode matrix includes control probes as extra rows.
# Split them off so QC can report on them without them polluting the
# gene expression matrix used for clustering.
separate_controls <- function(counts) {
  is_control <- grepl(paste(CONTROL_PREFIXES, collapse = "|"), rownames(counts))
  list(genes = counts[!is_control, , drop = FALSE],
       controls = counts[is_control, , drop = FALSE])
}

load_run <- function(h5_path, cells_csv_path) {
  counts <- Read10X_h5(h5_path)
  # Newer Seurat versions return one matrix per feature-type modality
  # (Gene Expression / Negative Control Probe / ...) instead of a single
  # combined matrix; recombine so separate_controls can split by rowname.
  if (is.list(counts)) counts <- do.call(rbind, counts)
  split <- separate_controls(counts)
  cells_meta <- read.csv(cells_csv_path)
  rownames(cells_meta) <- cells_meta$cell_id
  list(counts = split$genes, controls = split$controls, meta = cells_meta)
}

# One Seurat object per resolved sample, tagged with run_id/sample_id and
# carrying the spatial centroid + control-probe counts in metadata.
build_seurat <- function(counts, controls, meta, run_id, sample_id, cell_ids) {
  counts <- counts[, cell_ids, drop = FALSE]
  seu <- CreateSeuratObject(counts = counts, assay = "Xenium")
  seu <- AddMetaData(seu, meta[cell_ids, , drop = FALSE])
  seu$nCount_controls <- Matrix::colSums(controls[, cell_ids, drop = FALSE])
  seu$run_id <- run_id
  seu$sample_id <- sample_id
  seu
}

resolve_run <- function(row) {
  d <- load_run(row$h5_path, row$cells_csv_path)

  if (row$split_method == "single") {
    return(list(build_seurat(d$counts, d$controls, d$meta,
                              row$run_id, row$sample_id, colnames(d$counts))))
  }

  if (row$split_method == "annotation") {
    # coordinates.csv has a couple of "#" comment lines (selection names,
    # areas) before the real header; comment.char skips those.
    poly <- read.csv(row$annotation_csv_path, comment.char = "#")
    sel_names <- unique(poly$Selection)

    inside <- vapply(sel_names, function(s) {
      p <- poly[poly$Selection == s, ]
      sp::point.in.polygon(d$meta$x_centroid, d$meta$y_centroid, p$X, p$Y) > 0
    }, logical(nrow(d$meta)))

    n_hits <- rowSums(inside)
    if (any(n_hits > 1)) {
      warning(sum(n_hits > 1), " cells in ", row$annotation_csv_path,
              " fall inside more than one selection polygon; dropping them")
    }
    d$meta$sample_id <- NA_character_
    single_hit <- n_hits == 1
    d$meta$sample_id[single_hit] <- sel_names[max.col(inside[single_hit, , drop = FALSE])]
    d$meta <- d$meta[!is.na(d$meta$sample_id), ]

    groups <- split(rownames(d$meta), d$meta$sample_id)
    return(lapply(names(groups), function(g)
      build_seurat(d$counts, d$controls, d$meta, row$run_id, g, groups[[g]])))
  }

  if (row$split_method == "spatial_kmeans") {
    coords <- d$meta[, c("x_centroid", "y_centroid")]
    km <- kmeans(coords, centers = 2)
    d$meta$sample_id <- paste0(row$run_id, "_s", km$cluster)
    plot(coords$x_centroid, coords$y_centroid, col = km$cluster, pch = 16,
         cex = 0.2, asp = 1, main = paste("Spatial split -", row$run_id))
    groups <- split(rownames(d$meta), d$meta$sample_id)
    return(lapply(names(groups), function(g)
      build_seurat(d$counts, d$controls, d$meta, row$run_id, g, groups[[g]])))
  }

  stop("Unknown split_method: ", row$split_method)
}

manifest <- read.csv("manifest.csv", stringsAsFactors = FALSE)
samples <- unlist(lapply(seq_len(nrow(manifest)), function(i)
  resolve_run(manifest[i, ])), recursive = FALSE)
names(samples) <- vapply(samples, function(s) unique(s$sample_id), character(1))

cat("Resolved", length(samples), "samples from", nrow(manifest), "runs\n")
saveRDS(samples, "samples_list_raw.rds")


#read rds file
count_meta <- readRDS("samples_list_raw.rds")

samples <- readRDS("samples_list_raw.rds")
s <- samples[["S25-00422-1"]]

# metadata alone, as a data frame
md <- s@meta.data
View(md)                      # RStudio spreadsheet-style viewer

# metadata + specific genes' expression together
cell_table <- FetchData(s, vars = c("sample_id", "x_centroid", "y_centroid",
                                    "total_counts", "A2M", "ACE2"))
View(cell_table)

# export to CSV for Excel
write.csv(md, "S25-00422-1_meta.csv", row.names = TRUE)

