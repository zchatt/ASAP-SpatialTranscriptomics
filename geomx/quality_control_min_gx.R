# Scripts for reading in and performing quality control of GeoMx NGS data
# Modified from https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html. Thank you Nanostring!

# Libraries
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(knitr)
library(dplyr)
library(ggforce)
library(ggnewscale)
library(WGCNA)
library(ggplot2)
library(ggsignif)
library(readxl)
library(plyr)
library(scales)
library(gridExtra)
library(openxlsx)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(viridis)

# check if correct packages loaded
if(packageVersion("GeomxTools") < "2.1" & 
   packageVersion("GeoMxWorkflows") >= "1.0.1"){
  stop("GeomxTools and Workflow versions do not match. Please use the same version. 
    This workflow is meant to be used with most current version of packages. 
    If you are using an older version of Bioconductor please reinstall GeoMxWorkflows and use vignette(GeoMxWorkflows) instead")
}
if(packageVersion("GeomxTools") > "2.1" & 
   packageVersion("GeoMxWorkflows") <= "1.0.1"){
  stop("GeomxTools and Workflow versions do not match. Please use the same version, see install instructions above.")
}

############################################################################################
#### Inputs
############################################################################################
# # Locally produced data from human
run_name = "geomx_oct2023_min"
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min"
dcc_file_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/DCC_Files"
PKCFiles <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/pkcs/Hs_R_NGS_WTA_v1.0.pkc"
SampleAnnotationFile <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/Annotations_nm.xlsx" #Ensure annotations file correctly formatted; Tab = "Sample_ID" & sorted by Sample_ID

# # # Locally produced data from human
# run_name = "geomx_sep2023"
# analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis"
# dcc_file_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/DCC_Files"
# PKCFiles <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/pkcs/Hs_R_NGS_WTA_v1.0.pkc"
# SampleAnnotationFile <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/Annotations_remapped.xlsx"

# # Locally produced data from Substantia Nigra
# run_name = "CHA12467"
# analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/analysis"
# dcc_file_dir <- "//Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/DCC_Files"
# PKCFiles <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/pkcs/Hs_R_NGS_WTA_v1.0.pkc"
# SampleAnnotationFile <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/Annotations.xlsx"

# # # Nanotring produced data from Cortex
# run_name = "hu_brain_001"
# analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/analysis"
# dcc_file_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/workflow/dcc/run1_DCC"
# PKCFiles <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/workflow/pkc/Hs_R_NGS_WTA_v1.0.pkc"
# SampleAnnotationFile <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/Annotations_run1.xls"

## thresholds
target_genes_detected_samples.frac = 0 # formally set as 0.1
min_genes_detected.frac = 0.01 # formally set as 0.1

############################################################################################
#### Run
############################################################################################
setwd(analysis_dir)

# automatically list files in each directory for use
DCCFiles <- dir(file.path(dcc_file_dir), pattern = ".dcc$", full.names = TRUE)

# load data
# NOTE: the data columns of interest need to be present within the SampleAnnotationFile. QC plots will be made for each.
data_cols_interest <- c("roi","aoi","area","Brainbank_ID","AOI","Brainregion","Brainregion_2","PMD hs","Archive ys-2023","DV200","IHC-score","Age","Sex","Fam","Diagnosis_Broad","Brainbank","Diagnosis")

demoData <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = SampleAnnotationFile,
                         phenoDataSheet = "Template",
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = data_cols_interest,
                         experimentDataColNames = c("panel"))

# # NOTE: tip on accessing data
# sData(demoData) # sequencing data metrics for each ROI
# fData(demoData) # feature data for each probe
# exprs(demoData) # expression data for each feature x ROI

# # protocol suggests shifting 0 to 1 to enable in downstream transformations.
# # however we find 2.27% of counts == 1 (CHA12467), 5.85% of counts == 1 (hu_brain_001). Hence removal of < 2 should be performed below.
# table(exprs(demoData) == 1)/ length(exprs(demoData) == 1) * 100

# Shift counts to one
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)

# assign data for visium selection
demoData_vis <- demoData

# Select Segment QC
# Default QC cutoffs are commented in () adjacent to the respective parameters
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (75%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 10000,     # Maximum counts observed in NTC well (10000); Note this will be highly dependent on seq depth
       minNuclei = 20,         # Minimum # of nuclei estimated (20)
       minArea = 1000          # Minimum segment area (1000)
       )

demoData <- setSegmentQCFlags(demoData, qcCutoffs = QC_params)        

# Define Modules Used
pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

# Visualize Segment QC
col_by <- "segment"
# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

# Save QC plots
qc_hist1 <- QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)
ggsave("qc_1_hist.trimmed.png",qc_hist1, device = "png")
qc_hist2 <- QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)
ggsave("qc_2_hist.stitched.png",qc_hist2, device = "png")
qc_hist3 <- QC_histogram(sData(demoData), "Aligned (%)", col_by, 75)
ggsave("qc_3_hist.aligned.png",qc_hist3, device = "png")
qc_hist4 <- QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
ggsave("qc_4_hist.saturation.png",qc_hist4, device = "png")
dplot <- sData(demoData)
dplot$area <- as.numeric(dplot$area)
qc_hist5 <- QC_histogram(dplot, "area", col_by,1000, scale_trans = "log10")
ggsave("qc_5_hist.area.png",qc_hist5, device = "png")

# aggregated for publication
arrange <- ggarrange(plotlist=list(qc_hist1,qc_hist2,qc_hist3,qc_hist4), nrow=2, ncol=2, widths = c(2,2))
ggsave("Supp_Fig1.png",arrange, device = "png")


# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(demoData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}

ggsave("qc_5b_hist.geomean.png",plt, device = "png")

# detatch neg_geomean columns ahead of aggregateCounts call
pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(demoData)$NTC),
      col.names = c("NTC Count", "# of Segments"))

# plot all of the QC Summary information in a table.
qcss <- tableGrob(QC_Summary)
ggsave("qc_6_table.qcsummary.png",qcss,device = "png")

# Remove flagged segments
demoData <- demoData[, QCResults$QCStatus == "PASS"]
dim(demoData)

# Probe QC
# Set Probe QC Flags
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

# Exclude Outlier Probes
#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- subset(demoData, 
           fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

dim(ProbeQCPassed)
demoData <- ProbeQCPassed 

# Create Gene-level Count Data
# Check how many unique targets the object has
length(unique(featureData(demoData)[["TargetName"]]))

# collapse to targets
gxdat <- aggregateCounts(demoData)
dim(gxdat)

## Check if TYR present
fdat <- fData(demoData)
dat <- exprs(demoData)
gene_name <- "TYR"
gene <- dat[fdat$RTS_ID[fdat$SystematicName == gene_name],]

hist(gene/(sData(demoData)$NTC + 1), xlab = "TYR counts /(NTC counts + 1)", breaks = 50)


# Limit of Quantification
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(gxdat))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(gxdat)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(gxdat)[, vars[1]] * 
             pData(gxdat)[, vars[2]] ^ cutoff)
  }
}
pData(gxdat)$LOQ <- LOQ

## Filtering
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(gxdat)$Module == module
  Mat_i <- t(esApply(gxdat[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

table(LOQ_Mat["TYR",]) 

sData(demoData)[,"Aligned"]

x <- colSums(exprs(demoData), na.rm =T)[names(LOQ_Mat["TYR",colnames(tmp)])[LOQ_Mat["TYR",colnames(tmp)]]]
y <- colSums(exprs(demoData), na.rm =T)[names(LOQ_Mat["TYR",colnames(tmp)])[!LOQ_Mat["TYR",colnames(tmp)]]]
t.test(x,y)

# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(gxdat)$TargetName, ]

# Segment Gene Detection
# Save detection rate information to pheno data
pData(gxdat)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(gxdat)$GeneDetectionRate <-
  pData(gxdat)$GenesDetected / nrow(gxdat)

# Determine detection thresholds:<1%,2%,3%,4%, 5%, 10%, 15%, >15%
pData(gxdat)$DetectionThreshold <- 
  cut(pData(gxdat)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.02,0.03, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "2%","3%","4%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
qc_gd_bar <- ggplot(pData(gxdat),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

ggsave("qc_7_bar.genedetection.png",qc_gd_bar,device = "png")


# remove segments with less than min genes detected threshold
gxdat <- gxdat[, pData(gxdat)$GeneDetectionRate >= min_genes_detected.frac ]
dim(gxdat)

#Gene Detection Rate
# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(gxdat)]
fData(gxdat)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(gxdat)$DetectionRate <-
  fData(gxdat)$DetectedSegments / nrow(pData(gxdat))

# Evaluate genes of relevance
goi <- c("TYR")
goi_df <- data.frame(
  Gene = goi,
  Number_of_Segments = fData(gxdat)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(gxdat)[goi, "DetectionRate"]))

ss <- tableGrob(goi_df)
ggsave("qc_12_table.relevantgenes..png",ss,device = "png")

# Gene Filtering
# Plot detection rate:
plot_detect <- data.frame(Freq = c(0.5,1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.005,0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(gxdat)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(gxdat))
rownames(plot_detect) <- plot_detect$Freq

gdr_bar <- ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")
ggsave("qc_13_bar.genesdetected.png",gdr_bar,device = "png")

# Subset to target genes detected in at least X% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(gxdat), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
gxdat <- 
  gxdat[fData(gxdat)$DetectionRate >= target_genes_detected_samples.frac &
                    !(fData(gxdat)$TargetName %in% neg_probes), ]
dim(gxdat)

# Write samples that pass quality control to file
save(gxdat,file=paste0(run_name,"_gt",target_genes_detected_samples.frac,".gs",min_genes_detected.frac,"_qc.gx.Rdata"))


# #### selection of samples for visium analysis - done from n=901 data
# # scale DV200 and gene detected
# dplot <- sData(demoData)
# dplot$scaled_DV200<- scale(as.numeric(dplot$DV200))
# dplot$GenesDetected <- pData(gxdat)$GenesDetected
# dplot$scaled_GenesDetected <- scale(pData(gxdat)$GenesDetected)
# 
# # aggregate
# library(dplyr)
# mdat <- dplot %>%
#   group_by(Brainbank_ID) %>%
#   summarise_at(vars(scaled_DV200,scaled_GenesDetected), list(name = mean))
# 
# mdat <- as.data.frame(mdat)
# 
# var1 <- c("Brainbank_ID","Brainregion","PMD hs","Archive ys-2023", "Age", "Sex","Fam","Diagnosis Broad","Diagnosis","DV200","GenesDetected")
# tmp <- merge(mdat, dplot[,var1], by="Brainbank_ID", all.y = FALSE)
# tmp <- tmp[!duplicated(tmp$Brainbank_ID),]
# 
# tmp$dv200.gd_rank <- tmp$scaled_DV200_name + tmp$scaled_GenesDetected
# 
# tmp <- tmp[order(-tmp$dv200.gd_rank),]
# 
# write.table(tmp,file="summarised_geomx_qc_120923.txt",sep="\t",quote = F,row.names = F)
# 
# # manually select patients balancing for Dx (CTR, ePD, lPD and ILBD), Sex, brain region
# setwd(analysis_dir)
# tmp2 <- read_xlsx("summarised_geomx_qc_120923_select.xlsx",1)
# 
# tmp2[tmp2$visium_round1 == "Y" | tmp2$visium_round2 == "Y",] %>%
#   group_by(Diagnosis,Sex,Brainregion) %>%
#   dplyr::summarise(n = n())
# 
# tmp2[tmp2$visium_round1 == "Y" | tmp2$visium_round2 == "Y",] %>%
#   group_by(Diagnosis) %>%
#   dplyr::summarise(n = n())
# 
# tmp2[tmp2$visium_round1 == "Y" | tmp2$visium_round2 == "Y",] %>%
#   group_by(Brainregion) %>%
#   dplyr::summarise(n = n())
# 
# tmp2[tmp2$visium_round1 == "Y" | tmp2$visium_round2 == "Y",] %>%
#   group_by(Sex) %>%
#   dplyr::summarise(n = n())
# 
# tmp2[tmp2$visium_round1 == "Y" | tmp2$visium_round2 == "Y",] %>%
#   group_by(Brainbank_ID2,Diagnosis) %>%
#   dplyr::summarise(n = n())



##### tabulation of cohort statistics
# create a glm of cohort stats
library(stargazer)

dplot$Sex[dplot$Sex %in% c("male","Male")] <- "M"
dplot$Sex[dplot$Sex %in% c("female","Female")] <- "F"
dplot$Diagnosis <- as.factor(dplot$Diagnosis)
dplot$area <- as.numeric(dplot$area)
dplot$DV200 <- as.numeric(dplot$DV200)
dplot$Age <- as.numeric(dplot$Age)

lm.model <- glm(Diagnosis ~ Age + Sex + DV200 + area,family = binomial , data=dplot)

# print glm to file
stargazer(lm.model, title="Results", align=TRUE, out="star_descriptive.doc",type = "html")

## tabulate data for tables
tmp <- dplot[!duplicated(paste0(dplot$Brainregion,dplot$Brainbank)),]

# cases
tmp %>%
  group_by(Diagnosis,Brainregion) %>%
  dplyr::summarise(count = n_distinct(Brainbank))
tmp %>%
  group_by(Diagnosis) %>%
  dplyr::summarise(count = n_distinct(Brainbank))
dplot %>%
  group_by(Diagnosis) %>%
  dplyr::summarise(count = n_distinct(roi_id))
tmp %>%
  group_by(Diagnosis) %>%
  dplyr::summarise(count = n_distinct(Brainbank_ID_sub))

# Age
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(mean_age = mean(Age,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(sd_age = sd(Age,na.rm=TRUE))

#Sex
tmp %>%
  group_by(Diagnosis,Sex) %>%
  dplyr::summarise(count = n_distinct(Brainbank_ID))

# Dv200
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(mean_dv200 = mean(DV200,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(sd_dv200 = sd(DV200,na.rm=TRUE))

# Area
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(mean_area = mean(area,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(sd_dv200 = sd(DV200,na.rm=TRUE))


# ## tabulate data following ROI removal
# dplot <- pData(gxdat)
# col_names <- colnames(dplot)
# for(column in data_cols_interest) {
#   col_number = which(colnames(protocolData(gxdat)) == column)
#   dplot <- cbind(dplot,protocolData(gxdat)[[col_number]])
# }
# colnames(dplot) <- make.names(c(col_names,data_cols_interest))
# dplot$roi_id <- row.names(dplot)
# dplot$area <- as.numeric(dplot$area)
# dplot$DV200 <- as.numeric(dplot$DV200)
# dplot$Age <- as.numeric(dplot$Age)
# 
# # create DV200 and area bins
# dplot <- dplot %>% mutate(dv200_bin = cut(DV200, breaks=5))
# dplot <- dplot %>% mutate(area_bin = cut(area, breaks=5))
# 
# # cases
# dplot %>%
#   group_by(Diagnosis,Brainregion) %>%
#   dplyr::summarise(count = n_distinct(roi_id))