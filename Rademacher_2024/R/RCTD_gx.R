# RCTD analysis of GeoMx datasets for Rademacher et al., 2024


library(spacexr)
library(Seurat)
library(Matrix)
library(doParallel)
library(ggplot2)
library(plyr)
library(data.table)
library(msigdbr)
library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(EnhancedVolcano)
library(viridis)
library(rstatix)
library(ComplexHeatmap)
library(readxl)
library(lme4)
library(multcomp)
library(ggpubr)
source("/Users/zacc/github_repo/spacexr/R/CSIDE_plots.R")
register(SerialParam())


############################################################################################
#### Inputs
############################################################################################

analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
rdata = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata"
results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results_run180124'
contrast_path <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/contrasts_matrix.xlsx"
setwd(analysis_dir)

# single-cell data 
load("/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/merged_kam.sil.web_seurat.Rdata")

# NOTE; taking excitatory, inhibitory and non-neurons from Kamath, NE and 5-HT from Webber, DA neuron populations from Siletti re-code
toMatch <- c("Ex_","Inh_","Astro_","Endo_","Ependyma_","Macro_","MG_","Olig_","OPC_","5HT","NE",
             names(table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Siletti"])))
group1 <- grep(paste(toMatch,collapse="|"), merge.combined.sct@meta.data$cell_type_merge)
sc_obj <- merge.combined.sct[,group1]


############################################################################################
###### Part 1: Create ground-truths from snRNAseq datasets used for decon
############################################################################################
### create ground-truths from snRNAseq data

kamath_all <- table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Kamath"])
toMatch <- c("Ex_","Inh_","Astro_","Endo_","Ependyma_","Macro_","MG_","Olig_","OPC_")
group1 <- grep(paste(toMatch,collapse="|"), names(kamath_all))
kamath_nonDA <- kamath_all[group1]
kamath_nonDA_truth <- kamath_nonDA/ sum(kamath_nonDA)

kamath_all_truth <- table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Kamath"]) /
  sum(table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Kamath"]))


############################################################################################
###### Part 2:  Format geomx spatial data
############################################################################################
# setwd
setwd(results_folder)

# load geomx normalised data
load(rdata)
# select samples
#keep_index <- gxdat_s$Diagnosis == "CTR"
keep_index <- gxdat_s$Diagnosis != "NTC"
# load in counts matrix
#counts <- gxdat_s@assays$GeoMx@counts[,keep_index]  # load in counts matrix
counts <- gxdat_s@assays$RNA@counts[,keep_index]  # load in counts matrix
# create metadata matrix 
meta <- gxdat_s@meta.data[keep_index,]
# confirm names
table(colnames(counts) == row.names(meta))
# load in coordinates
coords = as.data.frame(gxdat_s@reductions$umap@cell.embeddings[keep_index,]) # we are using UMAP coordinates as dummy variables until we obtain DSP ROI locations
colnames(coords) <- c("imagerow","imagecol")
# process counts
nUMI <- colSums(counts) # In this case, total counts per spot


### ii) deconvolute using Kamath
set.seed(123)
# construct ST to deconvolute for SN tissue
group2 <- meta$ROI %in% c("SND","SNL","SNM","SNV","VTA","LC","RN")
puck <- SpatialRNA(coords[group2,], counts[,group2], nUMI[group2])
barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')

# construct references
# NOTE; taking all cells from Kamath
toMatch <- c(names(table(merge.combined.sct@meta.data$cell_type_merge[merge.combined.sct@meta.data$dataset_merge == "Kamath"])))
group1 <- merge.combined.sct@meta.data$cell_type_merge %in% toMatch
sc_obj <- merge.combined.sct[ ,group1]

# remove cell-types with < 25 cells AND/OR cells of interes
cells_keep <- names(table(sc_obj@meta.data$cell_type_merge)[table(sc_obj@meta.data$cell_type_merge) > 25])
#cells_keep <- c("CALB1_CALCR","CALB1_CRYM_CCDC68","CALB1_GEM", "CALB1_PPP1R17","CALB1_RBP4","CALB1_TRHR",
#"SOX6_AGTR1","SOX6_DDT","SOX6_GFRA2", "SOX6_PART1")
sc_obj <- sc_obj[,sc_obj@meta.data$cell_type_merge %in% cells_keep]
# get count data and restrict to only genes within GeoMx
counts_sc <- sc_obj[["RNA"]]$counts
counts_sc <- counts_sc[row.names(counts_sc) %in% row.names(puck@counts),]

# set cell types and quant nUMI
cell_types <- setNames(sc_obj@meta.data$cell_type_merge, colnames(sc_obj))
cell_types <- as.factor(cell_types)
nUMIsc <- colSums(counts_sc)
## create reference
reference <- Reference(counts_sc, cell_types, nUMIsc)

##  run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# add deconvolution to meta-data and save
table(row.names(normalize_weights(myRCTD@results$weights)) == row.names(meta))
meta_rctd <- cbind(meta,normalize_weights(myRCTD@results$weights))

############################################################################################
###### Part 3: Group-wise comparisons of DA neuron proportions
############################################################################################
# select just CTRs, ePD and lPD subjects
dat <- meta_rctd[meta_rctd$ROI == "SNV" & meta_rctd$segment == "Full ROI" & meta_rctd$Diagnosis %in% c("CTR", "ePD"),]

# sum DA neuron values
dat$da_neuron <- rowSums(dat[,c("SOX6_DDT","SOX6_GFRA2","SOX6_PART1", "CALB1_GEM",
                                "CALB1_PPP1R17","CALB1_RBP4","CALB1_TRHR","CALB1_CALCR",
                                "CALB1_CRYM_CCDC68")])

# calculate as % of controls
mean_controls <- median(dat$da_neuron[dat$Diagnosis == "CTR"])
dat$da_neuron_prc.ctr <- (dat$da_neuron / mean_controls) * 100

# format
data_table <- dat
colnames(data_table) <- make.names(colnames(data_table))
data_table <- data_table[,c("Diagnosis","da_neuron_prc.ctr","Brainbank_ID","Age","Sex","PMD.hs")]

# define variables
x_variable = "Diagnosis"
y_variable = "da_neuron_prc.ctr"
x_lab = "Diagnosis"
y_lab = "DA neurons (% of Controls)"
colour_palette = c("CTR"= "black","ePD" = "pink","lPD" = "red") 


# Fit a linear mixed-effects model
model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID) + (1 | Age) + (1 | Sex) + (1 | PMD.hs)")), data = data_table)

# perform post-hoc tests to compare different regions using Tukey method
posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
summary(posthoc)

# format for plotting
tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)

# make violin plot 
bxp <- ggviolin(
  data_table, x = x_variable, y = y_variable, 
  fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + ylab(y_lab) + xlab(x_lab) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                   geom="pointrange", color="black") +
  geom_hline(yintercept = 100, lty = 2) + scale_y_continuous(breaks=c(20,40,60,80,100,200,300)) +
  geom_hline(yintercept = c(20,40,60,80), lty = 2, color = "grey") 
  
ggsave(paste0("violin_ctr.epd_DA_rctd_fullroi",y_variable,".png"), bxp, width = 4, height = 4)
 

