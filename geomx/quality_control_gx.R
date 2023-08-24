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
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(readxl)
library(plyr)
library(scales)
library(gridExtra)

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
# # Locally produced data from Substantia Nigra
# run_name = "CHA12467"
# analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/analysis"
# dcc_file_dir <- "//Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/DCC_Files"
# PKCFiles <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/pkcs/Hs_R_NGS_WTA_v1.0.pkc"
# SampleAnnotationFile <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/Annotations.xlsx"

# # Nanotring produced data from Cortex
run_name = "hu_brain_001"
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/analysis"
dcc_file_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/workflow/dcc/run1_DCC"
PKCFiles <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/workflow/pkc/Hs_R_NGS_WTA_v1.0.pkc"
SampleAnnotationFile <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/Annotations_run1.xls"

## thresholds
target_genes_detected_samples.frac = 0.05
min_genes_detected.frac = 0.05

############################################################################################
#### Run
############################################################################################
setwd(analysis_dir)

# automatically list files in each directory for use
DCCFiles <- dir(file.path(dcc_file_dir), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)

# load data
# NOTE: the data columns of interest need to be present within the SampleAnnotationFile. QC plots will be made for each.
data_cols_interest <- c("aoi", "roi","Brainbank_ID","Diagnosis","Brainregion","DV200")
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

demoData <-setSegmentQCFlags(demoData, qcCutoffs = QC_params)        

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
qc_hist5 <- QC_histogram(sData(demoData), "area", col_by, 1000, scale_trans = "log10")
ggsave("qc_5_hist.area.png",qc_hist5, device = "png")


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
target_demoData <- aggregateCounts(demoData)
dim(target_demoData)

# Limit of Quantification
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_demoData)[, vars[1]] * 
             pData(target_demoData)[, vars[2]] ^ cutoff)
  }
}
pData(target_demoData)$LOQ <- LOQ

## Filtering
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_demoData)$Module == module
  Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]

# Segment Gene Detection
# Save detection rate information to pheno data
pData(target_demoData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
  pData(target_demoData)$GenesDetected / nrow(target_demoData)

# Determine detection thresholds:<1%,2%,3%,4%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <- 
  cut(pData(target_demoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.02,0.03, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "2%","3%","4%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
qc_gd_bar <- ggplot(pData(target_demoData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

ggsave("qc_7_bar.genedetection.png",qc_gd_bar,device = "png")

# plots prior to segment removal
dplot <- pData(target_demoData)
col_names <- colnames(dplot)
for(column in data_cols_interest) {
  col_number = which(colnames(protocolData(target_demoData)) == column)
  dplot <- cbind(dplot,protocolData(target_demoData)[[col_number]])
}
colnames(dplot) <- c(col_names,data_cols_interest)

a <- labels2colors(dplot$Brainregion)
a <- a[order(dplot$GeneDetectionRate)]
qc_gd_1 <- ggplot(dplot, aes(x = reorder(roi, GeneDetectionRate), y = GeneDetectionRate*100, fill = segment)) + 
  geom_bar(stat = "identity") + ylab("Gene Detection Rate (%)") + xlab("ROI") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = a, size=5)) + ggtitle("xlabels = Brainregion") +
  geom_hline(yintercept = min_genes_detected.frac * 100)
ggsave("qc_8_bar.genedetection.png",qc_gd_1,device = "png")

a <- labels2colors(dplot$Brainbank_ID)
a <- a[order(dplot$GeneDetectionRate)]
qc_gd_2 <-ggplot(dplot, aes(x = reorder(roi, GeneDetectionRate), y = GeneDetectionRate*100, fill = segment)) + 
  geom_bar(stat = "identity") + ylab("Gene Detection Rate (%)") + xlab("ROI") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = a, size=5)) + ggtitle("xlabels = Brainbank_ID") +
  geom_hline(yintercept = min_genes_detected.frac * 100)
ggsave("qc_9_bar.genedetection.png",qc_gd_2,device = "png")

a <- labels2colors(dplot$Diagnosis)
a <- a[order(dplot$GeneDetectionRate)]
qc_gd_3 <-ggplot(dplot, aes(x = reorder(roi, GeneDetectionRate), y = GeneDetectionRate*100, fill = segment)) + 
  geom_bar(stat = "identity") + ylab("Gene Detection Rate (%)") + xlab("ROI") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = a, size=5)) + ggtitle("xlabels = Diagnosis") +
  geom_hline(yintercept = min_genes_detected.frac * 100)
ggsave("qc_10_bar.genedetection.png",qc_gd_3,device = "png") 

qc_gd_4 <- ggplot(dplot, aes(x=as.character(DV200), y = GeneDetectionRate*100, fill=segment)) + 
  geom_boxplot() + ylab("Gene Detection Rate (%)") + 
  stat_compare_means( label = "p.format", method="wilcox.test", hide.ns=T, paired=F)
ggsave("qc_11_box.genedetection.png",qc_gd_4,device = "png")


# remove segments with less than 5% of the genes detected.
target_demoData <-
  target_demoData[, pData(target_demoData)$GeneDetectionRate >= min_genes_detected.frac ]

dim(target_demoData)

#Gene Detection Rate
# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
  fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))

# Test genes of relevance to Brain and PD
# Using Neuron markers DaN markers as well as top DEGs identified across cell-types in Smajić et al., Brain, 2022
goi <- c("CALB1","TH","CADPS2",
         "ST6GALNAC3", # "RP11-556E13.1","NEAT1", # non-coding RNA
         "MAN1C1", "S100A6","PLXND1", "FMO5", "SNX31","TFAP2D","HSPA1A","CALD1","SLC44A1", #"HIST1H2AC" == H2AC6
         "H2AC6","CFAP157", "SLC4A11", "CIITA")
goi_df <- data.frame(
  Gene = goi,
  Number_of_Segments = fData(target_demoData)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_demoData)[goi, "DetectionRate"]))

ss <- tableGrob(goi_df)
ggsave("qc_12_table.relevantgenes..png",ss,device = "png")

# Gene Filtering
# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData))
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
negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_demoData <- 
  target_demoData[fData(target_demoData)$DetectionRate >= target_genes_detected_samples.frac |
                    fData(target_demoData)$TargetName %in% neg_probes, ]
dim(target_demoData)

# Write samples that pass quality control to file
save(target_demoData,file=paste0(run_name,"_qc.gx.Rdata"))
