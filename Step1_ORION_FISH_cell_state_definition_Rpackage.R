# install.packages("remotes")
#remotes::install_github("clemenshug/cycif.importer")

library(cycif.importer) # this is the package CH wrote to import cell feature table from McMicro
#the package: https://github.com/clemenshug/cycif.importer 
library(dplyr)
library(tidyverse)
library(sfnetworks)



# Define double gating pairs (optional - defaults will be used if not specified)
double_gates_df <- data.frame(
  marker1 = c("CD8a",  "CD4", "CD11c", "CD11c",  "CD4", "CD8a", "CD68"),
  marker2 = c("Ki67", "Ki67", "CD68", "TIM3",  "TIM3", "TIM3", "TIM3"),
  stringsAsFactors = FALSE
)

marker_combo_immune_met_ORIONFISH_2026 <- list(
  
  
  # APC combinations and epithelial cell combinations
  
  CD11cpCD68p = ~CD11cp & CD68p,
  CD11cpCD68n = ~CD11cp & !CD68p,
  EcadpPAX8p = ~`E-cadherinp`, #epithelial cells
  PAX8pEcadherinpKi67p = ~((PAX8p & Ki67p) | (`E-cadherinp` & Ki67p)), #epithelial proliferating cells
  PAX8pEcadherinpLDHp = ~((PAX8p & LDHp) | (`E-cadherinp` & LDHp)), #epithelial cells with high LDH (metabolic marker)
  PAX8pEcadherinpgH2axp = ~((PAX8p & gH2axp) | (`E-cadherinp` & gH2axp)), #epithelial cells with DNA damage (gH2ax is a marker of DNA damage)
  
  # Activated T/DC cells (complex expression)
  
  CD8apPD1pKi67n = ~CD8ap & PD1p & !Ki67p,
  CD8apPD1pKi67p = ~((CD8ap & PD1p) | (CD8ap & Ki67p)), #general activation 
  CD4pPD1pKi67p = ~((CD4p & PD1p) | (CD4p & Ki67p)), #general activation 
  CD4pPD1pKi67n = ~CD4p & PD1p & !Ki67p
)

rois <- cycif_load_roi_data("ROI")

#check the variables under Name
unique(rois$Name)

ROI_resolved <- resolve_roi_overlaps(rois)




#Run the analysis pipeline: 

result <- cycif_pipeline(
  counts = "quantification_ORIONFISH", # quantification folder
  rois = ROI_resolved, 
  gate_table = "gateTable.csv", # file in the working directory
  slide_info = "slideINFOm.csv", # clinical metadata, file in the working directory
  roi_info = "tableROI.csv", #file in the working directory
  output_dir = "output_ORIONFISH",
  #roi_priority = c("p53signatures", "normal FT"), # Has to match the "Name" column in the ROI df #optional
  summary_filters = list(
    `PanCK+` = ~`E-cadherinp`, # here we are defining the epithelial marker or tumor marker
    `PanCK-` = ~!`E-cadherinp` # ditto- negative of that marker will serve as stroma
  ),
  summary_groups = list(
    #`by_brca` = c("BRCA", "ROIcat"),# if we want to group by BRCA mut and WT
    `by ROI_type`= c("slideName", "ROIcat") # columns directly contribute to the output file
    #`by ROI` = c("slideName", "ROI", "ROIcat"), # name from the annotated ROIs files from omero
    #`by slide` = c("slideName") #slideName needs to be consistent for matching ROIs, gate, cell feature table
  ),
  #slide_filter = c("LSP11064", "LSP11076"), #default is all samples but if we have to run few samples we can filter here
  scale_factor = 0.325,  # Scale factor for pixel to micron conversion #unbinned/binned
  sampling_mode = "roi_only",  # Sampling mode: "all_cells" or "roi_only"
  sample_size = NULL,  # Sample size for resampling, set NULL to include all cells
  expand_distance = 0,  # Minimum distance for ROI expansion (microns)
  double_gates = double_gates_df,
  marker_combinations = marker_combo_immune_met_ORIONFISH_2026,  # when have mix of populations, run this line #
  cell_file_pattern = "*-quantification-fish.csv*",
  roi_file_pattern = "*-rois.csv"
  #marker_combinations = c(
  #get_default_marker_combinations(),
  #list(
  #abc_cells = ~asdfp | asdfasp
)
#)
#)

#read zipped files: # or simply install 7-zip to extract files on the desktop
readsampledcells <- read_csv("sampled_cells.csv.gz")
ROI_type_panCKpos <- read_csv("summary_by ROI_type__panCK+.csv.gz")
ROI_type_panckneg <- read_csv("summary_by ROI_type__panCK-.csv.gz")



########## IMPORTANT

#If QC needed, using single cell file, 
#then follow the QC script first then generate panckpos and neg files from the QC final quant file




