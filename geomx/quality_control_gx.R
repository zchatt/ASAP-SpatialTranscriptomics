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
library(plyr )

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
  
  # to remove current package version
  # remove.packages("GeomxTools")
  # remove.packages("GeoMxWorkflows")
  # see install instructions above 
}

# inputs
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467"
dcc_file_dir <- "/Users/zacc/BaseSpace/CHA12467-ZC-6214210/GeoMx_NGS_Pipeline_08_01_2023_10_52_16-9375434/GeoMx_NGS_Pipeline_08_01_2023_10_52_16-ds.a1b6d6a5800c41b5a21056c3289f6943"
PKCFiles <- "/Users/zacc/BaseSpace/CHA12467-ZC-6214210/GeoMx_NGS_Pipeline_08_01_2023_10_52_16-9375434/GeoMx_NGS_Pipeline_08_01_2023_10_52_16-ds.a1b6d6a5800c41b5a21056c3289f6943/pkcs/Hs_R_NGS_WTA_v1.0.pkc"
SampleAnnotationFile <- "/Users/zacc/USyd/spatial_transcriptomics/GeoMx_methods/GeoMx_round_1_YH_20230719T2320/Annotations.xlsx"

#########
## RUN ##
#########
setwd(analysis_dir)

# automatically list files in each directory for use
DCCFiles <- dir(file.path(dcc_file_dir, "DCC_Files"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)

# load data
# NOTE: the data columns of interest need to be present within the SampleAnnotationFile. QC plots will be made for each.
data_cols_interest <- c("aoi", "roi","Brainbank_ID","Diagnosis","Brainregion")
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


# # within the protocol they suggest shifting 0 to 1 to enable in downstream transformations.
# # however we find 2.27% of counts == 1. Hence removal of < 2 should be performed below.
# table(exprs(demoData) == 1)/ length(exprs(demoData) == 1) * 100

# Shift counts to one
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)

################
## Segment QC ##
################

# Select Segment QC
# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (75%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 10000,     # Maximum counts observed in NTC well (10000); Note this will be highly dependent on seq depth
       #minNuclei = 20,         # Minimum # of nuclei estimated (20)
       #minArea = 1000          # Minimum segment area (1000)
       )
demoData <-
  setSegmentQCFlags(demoData, 
                    qcCutoffs = QC_params)        

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
QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)
QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)
QC_histogram(sData(demoData), "Aligned (%)", col_by, 75)
QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
QC_histogram(sData(demoData), "area", col_by, 1000, scale_trans = "log10")





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


#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(demoData, 
         fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
#> Features  Samples 
#>    18641      229
demoData <- ProbeQCPassed 



# Check how many unique targets the object has
length(unique(featureData(demoData)[["TargetName"]]))
#> [1] 18504

# collapse to targets
target_demoData <- aggregateCounts(demoData)
dim(target_demoData)
exprs(target_demoData)[1:5, 1:2]


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


# Save detection rate information to pheno data
pData(target_demoData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
  pData(target_demoData)$GenesDetected / nrow(target_demoData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <- 
  cut(pData(target_demoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

pData(target_demoData)$DetectionThreshold <- 
  cut(pData(target_demoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.02,0.03, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "2%","3%","4%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_demoData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")


## write gene detection rates to file
dplot <- pData(target_demoData)
dplot$roi <- protocolData(demoData)[["roi"]]
dplot$Brainbank_ID <- protocolData(demoData)[["Brainbank_ID"]]
dplot$Brainregion <- protocolData(demoData)[["Brainregion"]]
dplot$Diagnosis <- protocolData(demoData)[["Diagnosis"]]
#dplot <- dplot[pData(target_demoData)$GeneDetectionRate < 0.05,]
gcases_id <- read_excel("/Users/zacc/USyd/spatial_transcriptomics/GeoMx_methods/GeoMx_round_1_YH_20230719T2320/GeoMx_cases_zcformat.xlsx")
dplot$DV200 <- mapvalues(dplot$Brainbank_ID,gcases_id$Brainbank_ID,gcases_id$DV200)


pdf("geomx_genedetectionrate.pdf")
a <- labels2colors(dplot$Brainregion)
a <- a[order(dplot$GeneDetectionRate)]
ggplot(dplot, aes(x = reorder(roi, GeneDetectionRate), y = GeneDetectionRate*100, fill = segment)) + 
  geom_bar(stat = "identity") + ylab("Gene Detection Rate (%)") + xlab("ROI") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = a, size=5)) + ggtitle("xlabels = Brainregion")

a <- labels2colors(dplot$Brainbank_ID)
a <- a[order(dplot$GeneDetectionRate)]
ggplot(dplot, aes(x = reorder(roi, GeneDetectionRate), y = GeneDetectionRate*100, fill = segment)) + 
  geom_bar(stat = "identity") + ylab("Gene Detection Rate (%)") + xlab("ROI") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = a, size=5)) + ggtitle("xlabels = Brainbank_ID")

a <- labels2colors(dplot$Diagnosis)
a <- a[order(dplot$GeneDetectionRate)]
ggplot(dplot, aes(x = reorder(roi, GeneDetectionRate), y = GeneDetectionRate*100, fill = segment)) + 
  geom_bar(stat = "identity") + ylab("Gene Detection Rate (%)") + xlab("ROI") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = a, size=5)) + ggtitle("xlabels = Diagnosis")

dev.off()

pdf("boxplot_geomx_genedetectionrate.pdf")
ggplot(dplot, aes(x=DV200, y = GeneDetectionRate*100, fill=segment)) + 
  geom_boxplot() + ylab("Gene Detection Rate (%)") + 
  stat_compare_means( label = "p.format", method="wilcox.test", hide.ns=T, paired=F)
dev.off()




