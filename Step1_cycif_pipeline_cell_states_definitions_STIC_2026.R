# install.packages("remotes")
#remotes::install_github("clemenshug/cycif.importer")

library(cycif.importer) # this is the package CH wrote to import cell feature table from McMicro
library(dplyr)
library(tidyverse)


# Define double gating pairs (optional - defaults will be used if not specified)
double_gates_df <- data.frame(
  marker1 = c("CD163",  "CD68", "CD8a", "CD4", "CD4", "CD4", "CD4", "CD68", "CD8a", "CD8a", "PanCK", "PanCK", "HLA-E", "CD11c", "pTBK1", "CD11c", "Mx1", "CD11c", "PKM2", "PKM2", "ENO1", "ENO1", "LDH", "LDH", "CD11c", "ORF1p", "ORF1p", "MCT1"),
  marker2 = c("CD11c", "CD11c", "CD103","CD68", "Tbet", "FoxP3", "HLA-DR", "HLA-DR", "TIM3", "Tbet", "gH2Ax", "Ki67", "Mx1", "HLA-DR", "HLA-E", "TIM3", "pTBK1", "CD40", "LDH", "ENO1", "LDH", "FASN", "FASN", "Ki67", "MCT4", "Mx1", "HLA-E","LDH"),
  stringsAsFactors = FALSE
)

marker_combinations_exp_stic_2026 <- list(
  # Basic T cell combinations
  CD8apCD103n = ~CD8ap & !CD103p, #CTL
  CD4pPD1pKi67p = ~CD4p & (`PD1-repeatp` | Ki67p), #activated CD4
  CD4pTIM3p = ~CD4p & TIM3p, #terminally exh CD4
  CD4pTIM3pLAG3p = ~CD4p & (TIM3p | LAG3p), #CD4 dysfunction
  CD4pLAG3pFoxP3n = ~CD4p & LAG3p & !FoxP3p,   #Type1 regulatory CD4+ T cells (not Tregs)
  CD8apLAG3pTIM3n = ~CD8ap & LAG3p & !TIM3p, #only CD8+LAG3+
  CD8apLAG3pTIM3p = ~((CD8ap & `LAG3p`) | (CD8ap & TIM3p)), #total exhausted CD8+ T
  
  
  #NKG2A/Tbet T cell subtype
  NKG2ApCD8anCD4n = ~NKG2Ap & !CD8ap & !CD4p, # NK cell
  NKG2ApCD8apCD103n = ~NKG2Ap & CD8ap & !CD103p, # NKG2A+ CTL
  NKG2ApCD8apCD103p = ~NKG2Ap & CD8ap & CD103p, #NKG2A+ Trm
  CD8apTbetpNKG2ApCD103p = ~CD8ap & Tbetp & NKG2Ap & CD103p, #Tbet+NKG2A+Trm
  CD8apTbetpNKG2ApCD103n = ~CD8ap & Tbetp & NKG2Ap & !CD103p, #Tbet+NKG2A+ CTL
  
  
  # IFN pathway (OR of multiple markers)
  IFNp = ~pTBK1p | Mx1p |`HLA-Ep`,
  
  # APC combinations
  
  CD11cpCD68nCD103p = ~CD11cp & CD103p & !CD68p, #cDC1
  CD11cpCD68pCD103n = ~CD11cp & CD68p & !CD103p, #macrophage derived APCs
  CD11cpCD68nCD103n = ~CD11cp & !CD68p & !CD103p, #other dendritic cells
  CD11cpCD103nHLADRpCD68p = ~CD11cp & !CD103p & (`HLA-DRp` | CD40p) & CD68p, #macrophage derived activated dendritic cells
  CD11cpCD103pHLADRpCD40pCD68n = ~CD11cp & CD103p & (`HLA-DRp` | CD40p) & !CD68p, #activated cDC1
  CD11cpCD103nHLADRpCD40pCD68n = ~CD11cp & !CD103p & (`HLA-DRp` | CD40p) & !CD68p,  #activated dendritic cells
  CD11cpCD68nCD103pTIM3p = ~CD11cp & CD103p & !CD68p & TIM3p, #TIM3+ cDC1
  CD11cpCD68pCD103nTIM3p = ~CD11cp & CD68p & !CD103p & TIM3p, #TIM3+ macrophage derived APCs
  CD11cpCD68nCD103nTIM3p = ~CD11cp & !CD68p & !CD103p & TIM3p, #TIM3+other dendritic cells
  
  # Activated T/cDc1 cells (complex expression)
  CD8apPD1pKi67pLAG3nTIM3n = ~((CD8ap & `PD1-repeatp`) | (CD8ap & Ki67p)) & !LAG3p & !TIM3p,
  CD8apTbetpNKG2ApCD103pPD1pKi67pLAG3nTIM3n = ~CD8ap & Tbetp & NKG2Ap & CD103p & (`PD1-repeatp` | Ki67p) & !LAG3p & !TIM3p, # activated Tbet+NKG2A+ trm
  CD8apTbetpNKG2ApCD103nPD1pKi67pLAG3nTIM3n = ~CD8ap & Tbetp & NKG2Ap & !CD103p & (`PD1-repeatp`| Ki67p) & !LAG3p & !TIM3p, # activated Tbet+NKG2A+ CTL
  CD8apTbetpNKG2ApCD103pLAG3pTIM3p = ~CD8ap & Tbetp & NKG2Ap & CD103p & (LAG3p | TIM3p), #exh Tbet+CD8+ Trm
  CD8apTbetpNKG2ApCD103nLAG3pTIM3p = ~CD8ap & Tbetp & NKG2Ap & !CD103p & (LAG3p | TIM3p), #exh Tbet+CD8+ CTL
  
  #acidosis niche
  Lactp = ~LDHp & (ENO1p | PKM2p), # includes LDH and either ENO1 or PKM2, which are key glycolysis enzymes
  Glyp = ~PKM2p | ENO1p # any glycolysis enzyme #last line does not need this `,`
  
)


rois <- cycif_load_roi_data("ROI")

ROI_resolved <- resolve_roi_overlaps(rois)

#####if I check the slide, I need to make sure I change the df name#### Or ROI_resolved wont contain all slides######
#if I have to check any particular slide: 
#ROI_resolved_2 <- resolve_roi_overlaps(rois %>% filter(slideName == "LSP35353_P110"))

#plotting the ROIs including resolved: 

#plot_roi_sf(roi_df_to_sf(ROI_resolved) %>% filter(slideName == "LSP35328_P110"))



#Run the analysis pipeline: 

result <- cycif_pipeline(
  counts = "quantification", # quantification folder
  rois = ROI_resolved, 
  gate_table = "gateTable.csv", # file in the working directory
  slide_info = "slideINFOm.csv", # clinical metadata, file in the working directory
  roi_info = "tableROI.csv", #file in the working directory
  output_dir = "output",
  #roi_priority = c("p53signatures", "normal FT"), # Has to match the "Name" column in the ROI df #optional
  summary_filters = list(
    #epithelial = ~p53p | PanCKp,
    `PanCK+` = ~`PanCKp`, # here we are defining the epithelial marker or tumor marker
    `PanCK-` = ~!`PanCKp` # ditto- negative of that marker will serve as stroma
  ),
  summary_groups = list(
    #`by_brca` = c("BRCA", "ROIcat"),# if we want to group by BRCA mut and WT
    #`by ROI_type`= c("slideName", "ROIname") # columns directly contribute to the output file #what we can usually do
    #`by ROI` = c("slideName", "ROI", "ROIcat"), # name from the annotated ROIs files from omero
    `by slide` = c("slideName") #slideName needs to be consistent for matching ROIs, gate, cell feature table
  ),
  #slide_filter = c("LSP11064", "LSP11076"), #default is all samples but if we have to run few samples we can filter here
  scale_factor = 0.65,  # Scale factor for pixel to micron conversion #unbinned/binned 
  sampling_mode = "roi_only",  # Sampling mode: "all_cells" or "roi_only"
  sample_size = NULL,  # Sample size for resampling, set NULL to include all cells
  expand_distance = 0,  # Minimum distance for ROI expansion (microns)
  double_gates = double_gates_df,
  marker_combinations = marker_combinations_exp_stic_2026,  # when have mix of populations, run this line #
  cell_file_pattern = "*.csv*",
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
ROI_type_panCKpos <- read_csv("summary_by slide__PanCK+.csv.gz")
ROI_type_panckneg <- read_csv("summary_by slide__PanCK-.csv.gz")









