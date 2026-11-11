# stic-ms-2026
Gray Precancer Atlas v 2.0_KADER_Santagata Lab: important codes or new codes: 
Step1_ ORION-FISH or CyCIF: using cyccif.importer package to integrate gate table, ROIs, and quantifications from each tissue imaging platform. Also contains the definition of cell types and states.
CD8_subset_denominator_CD8_glmm: An example of calculating immune cell subset out of a particular type of immune cells compared to out of total cells. glmm stat test is incorporated. CD11c, CD4 subset can be caluclated the same way with small changes of immune populations. 

Xenium pipeline from raw data are run in order: scripts are numbered. 
`01_load_and_split.R` – load Xenium runs, split multi-specimen slides into samples
`02_qc.R` – per-sample QC filtering
 `normalisation/03_merge_normalize.R` – merge samples, PFlog1pPF normalisation
 `04_cluster_pflog.R` – PCA, Harmony, clustering on subset per sample
 `05_project_clusters_full.R` – transfer cluster labels to the full dataset
 `06_full_dataset_pca_umap.R` – PCA, Harmony, UMAP on the full dataset
 `11_histology_annotation_full.R` – assign lesion histology labels from annotation polygons
 This normalisation was only used for clustering, PCA, UMAP. 
 Neighborhood and Immune scoring was followed by Seurat normalisation. 
