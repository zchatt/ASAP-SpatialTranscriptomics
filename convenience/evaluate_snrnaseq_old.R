# scripts for 

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
library(lme4)
library(lme4)
library(multcomp)
library(randomForest)
library(e1071)
library(caret)
library(nnet)
library(edgeR)
library(Seurat)
library(data.table)
library(harmony)
library(RColorBrewer)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(enrichR)
library(writexl)


############################################################################################
### input
############################################################################################ 

setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min")

############################################################################################
### colour palettes
############################################################################################ 

Diagnosis_col = c("CTR"= "grey","ILBD" = "#00AFBB", "ePD" = "#E7B800","lPD" = "red")
age_group_col = magma(4)
ROI_col = c("SNL" = "darkorchid1",
            "SNV" = "purple",
            "SND" = "purple3",
            "SNM" = "purple4",
            "VTA" = "forestgreen",
            "LC" = "yellow",
            "SNV_young" = "pink")

Cell_col = c("SOX6_AGTR1" = brewer.pal(n = 9, name = "YlGnBu")[7] ,
            "SOX6_DDT" = brewer.pal(n = 9, name = "YlGnBu")[6] ,
            "SOX6_PART1" = brewer.pal(n = 9, name = "YlGnBu")[5],
            "SOX6_GFRA2" = brewer.pal(n = 9, name = "YlGnBu")[4],
            "CALB1_CALCR" = brewer.pal(n = 9, name = "YlOrRd")[2],
            "CALB1_TRHR" = brewer.pal(n = 9, name = "YlOrRd")[3],
            "CALB1_PPP1R17" = brewer.pal(n = 9, name = "YlOrRd")[4],
            "CALB1_CRYM_CCDC68" = brewer.pal(n = 9, name = "YlOrRd")[5],
            "CALB1_GEM" = brewer.pal(n = 9, name = "YlOrRd")[6] ,
            "CALB1_RBP4" = brewer.pal(n = 9, name = "YlOrRd")[7],
            "NE" = "yellow")

Cell_col2 = c("SOX6_AGTR1_TYRpos" = "red",
              "SOX6_AGTR1_TYRneg" = "purple",
              "NE" = "yellow")
              

############################################################################################
### Neuromelanin related gene sets
############################################################################################

skin_melanin_enzymes = c("TYR", "TYRP1","DCT") # note DCT = "TYRP2"
alt_PD = c("ALDH1A1", "ALDH3A1", "DDT", "CRABP1", "MAGED2","RBP4")
CA_precursor_NM_stock_trans_metab_A9.A10 = c("TH","DDC","COMT","ALDH","MAO","AR","ADH")
CA_precursor_NM_stock_trans_metab_A6 = c("TH", "DDC", "DBH", "MAO", "COMT", "MHPG")
CA_functional_DA = c("SLC18A2", "SLC6A3","SLC18A1")
CA_functional_NE = c("SLC18A2", "SLC18A1", "SLC6A2")
Stress_granules.free_radical_scavenging = toupper(c("DDX6","DDX1", "DDX3a", "DDX17","DDX3X", "PABPC1", "eIF3A", "eIF3B", "eIF4G1",
                                            "G3BP1","G3BP2", "RAB33A", "MAP1L3CB2","MAP1LC3B","NDUFS7", "ACOT8", "COA7", "PRDX2",
                                            "PRDX1","PRDX2","PRDX5","TUBA1B","GPN1","GPNMB","EGFR","BLVTB","BLVRB","GSTT1", "SIRT5"))
Lysosome_pathways = c("CTSD", "LAMP2", "MAP1LC3B", "SCARB2", "UBA52", "GBA1","GBA","FBXO16", "ATG5", "HSAPA4L", "HSAPA4","HSPA9", "UCHL1")

#enzymes upstream of euNM, DAQ-DAC-DHI (enzyme involved in these two steps and products, eg DDT)
upstream_euNM = c("TH","DDC","DBH","DDT")

# Genes linking skin pigmentation and Parkinson’s disease
genes_link_skinpig.PD = c("GCH1", "GPNMB", "HERC2", "LRRK2","MC1R","PRKN","SNCA", "TPCN2","TYR", "TRPM7", "VPS35")

# yuhong list genes interest
yf_genes <- c("DDT",
"MC1R","TYR","TYRP1",
"OCA2","ALDH1A1","ALDH1A3","CRABP1","MAGED2","RBP4",
"TH","DDC","DBH","MAO","MAOA","COMT","MHPG","VMAT2","DAT","SLC6A3","VMAT1","DDX6","DDX1","DDX3a",
"DDX17","PABPC1","elF3A","elF3B","elF4G1","G3BP1","G3BP2","RAB33A","MAP1L3CB2","NDUFS7",
"ACOT8","COA7","PRDX2","PRDX1","PRDX5","TUBA1B","GPN1","GPNMB","EGFR","BLVTB","GSTT1",
"SIRT5","CTSD","LAMP2","MAP1LC3B","SCARB2","UBA52","GBA1","FOXO16","FOXO1","ATG5","HSPA4L",
"HSAPA4","HSPA4","HSPA9","UCHL1")

# pigmentation network
pigmentation_network <- c("ADAM17","ADAMTS20","BRCA1","ED1","EDA","EDN3","EDNRB","EGFR",
"FGFR2","IKBKG","KIT","KITLG","KRT2A","KRT2","LMX1A","MCOLN3","MITF","PAX3","SFXN1","SNAI2",
"SOX10","SOX18","WNT1","WNT3A","DCT","GPNMB","MATP","SLC45A2","SLC45A2","RAB38","SILV","PMEL","TYR","TYRP1",
"AP3B1","AP3D1","VPS33A","CNO","BLOC1S4","HPS1","HPS3","HPS4","HPS5","HPS6","LYST","MU","GSTM1",
"OA1","GPR143","P","PLDN","BLOC1S6","RABGGTA","MLPH","MYO5A","MYO7A","RAB27A","ASIP","ATRN","GGT loci (several)","GGT1","GL","MC1R",
"MGRN1","POMC","ATP7A","ATP7B","BCL2","ERCC2","DCX","GSR","ITG2B","ITGA2B",
"ITGB3","MAP3K14","PH","PPY","C4A","C4B","U2AF1","ZFP362","ZNF362","MECP2","PITX2","RS1","SCO2","TYMP","MLC","MLC1")

# yuhong + pigmentation genes
#yf_pigment_genes <- unlist(read.delim(file="/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min/pigmentation+YuHong.genes.txt", header = F)[,1])
yf_pigment_genes <- unlist(read.delim(file="/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min/pigmentation+YuHong.genes_fixed.txt", header = F)[,1])
yf_pigment_genes <-  toupper(yf_pigment_genes)
yf_pigment_genes <- unique(yf_pigment_genes )

# other undefines
other <- c("ABCB6","ANKRD27","GLS","LAMP1","RAB32","RAB9A","MYEF2","SGSM2","FOXO1")


############################################################################################
### Evaluate DA neurons from Kamath et al.
############################################################################################

#----------- 1) get data -----------
public_dataset_location = "/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/" # this is the mapped project drive
working_directory = "/Users/zacc/USyd/spatial_transcriptomics/analysis/public_datasets"
setwd(working_directory)
cell_ranger_data = paste0(public_dataset_location,"kamath_2023/") # single-cell cell ranger output
# i) load single-cell data
sc_obj <- Read10X(cell_ranger_data) # single-cell cell ranger output
# ii) create seurat object
seurat_obj <- CreateSeuratObject(sc_obj)
rm(sc_obj)
# iii) update metadata
tsv <- read.delim(paste0(cell_ranger_data,"METADATA_PD.tsv")) # metadata
tsv <- tsv[-1,]
row.names(tsv) <- tsv$NAME
seurat_obj <- AddMetaData(object = seurat_obj, metadata = tsv)
umap_file <- list.files(cell_ranger_data, pattern="*UMAP.tsv", full.names=TRUE) # read in umap of cell IDs
cell_id <- do.call(rbind, lapply(umap_file,fread))
cell_id <- cell_id[!cell_id$NAME == "TYPE",]
row.names(cell_id) <- cell_id$NAME
seurat_obj <- AddMetaData(object = seurat_obj, metadata = cell_id)

# genes of interest
genes_to_plot <- yf_pigment_genes 

# get expression data
exp_dat <- as.matrix(seurat_obj[['RNA']]$counts[row.names(seurat_obj[['RNA']]$counts) %in% genes_to_plot,])
meta_dat_total <- as.data.frame(seurat_obj@meta.data)
#rm(seurat_obj)

# genes of interest
table(unique(genes_to_plot) %in% row.names(exp_dat))
genes_to_plot[!genes_to_plot %in% row.names(exp_dat) ]
genes_to_plot <- unique(genes_to_plot[genes_to_plot %in% row.names(exp_dat)])


#-----------  2) distribution of counts per gene highlighting TYR ----------- 

# get expression data and cell metadata
meta_dat <- meta_dat_total[meta_dat_total$Cell_Type %in% unique(meta_dat_total$Cell_Type)[grep("SOX6_|CALB1_",unique(meta_dat_total$Cell_Type))] &
                             meta_dat_total$Status == "Ctrl",]
exp_mat1 <- exp_dat[,row.names(meta_dat)]

# Summarize the total counts for each gene
gene_summaries <- exp_mat1 %>%
  as.data.frame() %>%  # Convert the matrix to a data frame
  rownames_to_column(var = "Gene") %>%  # Move gene names to a column
  mutate(total_count = rowSums(.[-1]))  # Calculate total counts for each gene

# Name of the gene you want to highlight
gene_of_interest <- "TYR"

# Extract the total count for the gene of interest
tyr_count <- gene_summaries %>%
  filter(Gene == gene_of_interest) %>%
  pull(total_count)

# Create a histogram using ggplot2
g1 <- ggplot(gene_summaries, aes(x = total_count + 1)) +
  #geom_histogram(binwidth = 0.1, fill = "grey", color = "black", alpha = 0.7) +
  geom_density(aes(y=..density..)) +
  geom_segment(aes(x = tyr_count, xend = tyr_count, y = 0, yend = 0.3), linetype = "dashed", color = "grey20") +
  labs(title = "Counts per gene Ctr SNpc DA neurons",
       x = "Log10(Total Counts + 1)",
       y = "Frequency") +
  scale_x_log10() + theme_bw() + 
  geom_label(label=gene_of_interest, x=log10(tyr_count),
    y=0.3, label.padding = unit(0.55, "lines"),  label.size = 0.35,color = "grey20",
    fill="white")

## Create boxplot TYR expression by cell-type
# Create a contingency table and Perform Fisher's Exact Test
x <- meta_dat$Cell_Type == "SOX6_AGTR1"
y <- exp_mat1["TYR",] >= 1
contingencyTable <- table(x,y)
testResult <- fisher.test(contingencyTable)
testResult 

# Calculate the percentage of each cell type with counts > 0 for TYR gene
gene_expression_df <- as.data.frame(exp_mat1)

cell_types_vector <- meta_dat$Cell_Type
tmp <- gene_expression_df %>%
  filter(row.names(exp_mat1) == "TYR") %>%
  pivot_longer(cols = everything(), names_to = "Cell", values_to = "Count") %>%
  mutate(HasCount = Count > 0) %>%
  group_by(CellType = cell_types_vector)

tmp2 <- table(tmp$HasCount,tmp$CellType)
dplot <- as.data.frame(cbind(Percentage = as.numeric(tmp2["TRUE",]/colSums(tmp2) * 100), CellType = colnames(tmp2)))
dplot$Percentage <- as.numeric(dplot$Percentage)
dplot$CellType <- factor(dplot$CellType, levels = dplot$CellType[order(-dplot$Percentage)])

g2 <- ggplot(dplot,aes( x = CellType, y = Percentage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percentage DA neurons > 0 TYR counts",
       x = "Cell Type",
       y = "% Cells > 0 TYR counts",
       fill = "Count > 0") +
  theme_minimal() + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  geom_text(label="***", x=1,y=1.8) + ylim(0,2)

arrange <- ggarrange(plotlist=list(g1,g2), nrow=2, ncol=2)
ggsave("plots_kamath.ctr.counts.pdf",arrange)

# check sequencing depth between cell-types
x = meta_dat$nCount_RNA[exp_dat["TYR",] >= 1]
y = meta_dat$nCount_RNA[exp_dat["TYR",] == 0]
t.test(x,y)

x = meta_dat$nCount_RNA[exp_mat1["TYR",] >= 1 & meta_dat$Cell_Type == "SOX6_AGTR1"]
y = meta_dat$nCount_RNA[exp_mat1["TYR",] == 0 & meta_dat$Cell_Type == "SOX6_AGTR1"]
t.test(x,y)

x = meta_dat$nCount_RNA[ meta_dat$Cell_Type != "SOX6_AGTR1"]
y = meta_dat$nCount_RNA[ meta_dat$Cell_Type == "SOX6_AGTR1"]
t.test(x,y)

x = meta_dat$nCount_RNA[ meta_dat$Cell_Type == "SOX6_AGTR1"]
y = meta_dat$nCount_RNA[ meta_dat$Cell_Type == "CALB1_CALCR"]
t.test(x,y)


#-----------  3) perform DEG analysis between each cell-type in CTR ----------- 

# format variables
meta_dat$Cell_Type <- as.character(meta_dat$Cell_Type)

# define cell-types
cell_types <- unique(meta_dat$Cell_Type)

# format covariates
meta_dat$Cell_Type <- as.factor(meta_dat$Cell_Type)
meta_dat$donor_id <- as.factor(meta_dat$donor_id)
meta_dat$Donor_PMI<- as.factor(meta_dat$Donor_PMI)
meta_dat$libname <- as.factor(meta_dat$libname)
meta_dat$sex <- as.factor(meta_dat$sex)
meta_dat$nCount_RNA <- as.numeric(meta_dat$nCount_RNA)
meta_dat$Donor_Age <- as.numeric(meta_dat$Donor_Age)

meta_dat$Broad_Cell_Type <- as.factor(unlist(lapply(strsplit(as.character(meta_dat$Cell_Type),"_"),function(x) x[1])))
broad_cell_types <- unique(meta_dat$Broad_Cell_Type)

## Run Voom
design <- model.matrix( ~ Cell_Type + donor_id + libname + Donor_Age + sex + nCount_RNA, data=meta_dat,droplevels = F)
colnames(design) <- make.names(colnames(design))

# fit design
v <- voom(exp_mat1,design)
vfit <- lmFit(v)

# get results and assign to list
res <- list()
cont.matrix <- makeContrasts(A="(Cell_TypeCALB1_CRYM_CCDC68 + Cell_TypeCALB1_GEM + Cell_TypeCALB1_PPP1R17 +          
                             Cell_TypeCALB1_RBP4 + Cell_TypeCALB1_TRHR + Cell_TypeSOX6_AGTR1 + Cell_TypeSOX6_DDT +                 
                             Cell_TypeSOX6_GFRA2 + Cell_TypeSOX6_PART1)",levels=design)
fit2 <- contrasts.fit(vfit, cont.matrix)
vfit2 <- eBayes(fit2)
options(digits=3)

# Select significant upregulated DEGs and assign to list
tmp <- topTable(vfit2,number = Inf,sort.by="P")
res[[1]] <- tmp[tmp$adj.P.Val < 0.05,]

# run through non-baseline covariates
for (val in 2:length(cell_types)){
  # Perform LIMMA contrasts
  cont.matrix <- makeContrasts(A=colnames(design)[val],levels=design)
  fit2 <- contrasts.fit(vfit, cont.matrix)
  vfit2 <- eBayes(fit2)
  options(digits=3)
  
  # Select significant upregulated DEGs and assign to list
  tmp <- topTable(vfit2,number = Inf,sort.by="P")
  res[[val]] <- tmp[tmp$adj.P.Val < 0.05 ,]
}

# name results
names(res) <- levels(meta_dat$Cell_Type)
res <- lapply(res, function(x) {
  x$Gene <- row.names(x)
  row.names(x) <- NULL
  return(x)
})

# write to file
write_xlsx(res, path ="Kamath_ctr_DEG_pigment.n130.xlsx")

# define any DEGs that are unique to a cell-type
check <- 1:length(res)
res_uniq <- list()
for (i in 1:length(check)){
  check_min <- check[check != i]
  tmp <- do.call(rbind,res[check_min])
  tmp2 <- as.data.frame(res[[i]])
  
  tmpa <- tmp[tmp$logFC > 0,]
  tmp2a <- tmp2[tmp2$logFC > 0,]
  tmp3 <- tmp2a[!tmp2a$Gene %in% tmpa$Gene,]
  
  tmpa <- tmp[tmp$logFC < 0,]
  tmp2a <- tmp2[tmp2$logFC < 0,]
  tmp4 <- tmp2a[!tmp2a$Gene %in% tmpa$Gene,]
  
  res_uniq[[i]] <- rbind(tmp3,tmp4)
}
names(res_uniq) <- levels(meta_dat$Cell_Type)
res_uniq <- do.call(rbind,res_uniq)
res_uniq$Cell_type <-  unlist(lapply(strsplit(row.names(res_uniq),"\\."), function(x) x[1]))

# calculate percentile of each gene
total_expression <- rowSums(exp_dat)
percentile_exp <- ecdf(total_expression)(total_expression) * 100
names(percentile_exp) <- row.names(exp_dat)
sort(percentile_exp[names(percentile_exp) %in% genes_to_plot ])

# calculate the number of cells detected in
num_cells_expressed <- apply(exp_dat, 1, function(x) sum(x > 0))
names(num_cells_expressed) <- row.names(exp_dat)
sort(num_cells_expressed[names(num_cells_expressed) %in% genes_to_plot ])

# merge
meta1 <- meta_dat
meta1 <- bind_cols(meta1, as.data.frame(t(exp_mat1)))

# select cols interest
meta1 <- meta1[,c("Cell_Type",genes_to_plot)]

# long format
meta <- pivot_longer(meta1, -Cell_Type, names_to="Gene", values_to="Expression")

# create summary
meta_summary <- meta %>%
  dplyr::group_by(Cell_Type, Gene) %>%
  dplyr::summarise(Avg = mean(Expression),
                   Pct = sum(Expression > 0) / length(Expression) * 100)

# create gene groups
meta_summary$group <- as.character(meta_summary$Gene)
meta_summary$group[meta_summary$Gene %in% yf_genes] <- "Genes interest"
meta_summary$group[meta_summary$Gene %in% other] <- "Genes interest"
meta_summary$group[meta_summary$Gene %in% pigmentation_network] <- "Pigmentation Network"
meta_summary$group[meta_summary$Gene %in% upstream_euNM] <- "Upstream euNM"
meta_summary$group[meta_summary$Gene %in% c(CA_precursor_NM_stock_trans_metab_A9.A10,CA_precursor_NM_stock_trans_metab_A6)] <- "CA precursor"
meta_summary$group[meta_summary$Gene %in% c(CA_functional_DA,CA_functional_NE)] <- "CA Functional"
meta_summary$group[meta_summary$Gene %in% Stress_granules.free_radical_scavenging] <- "Stress Granule/ Scavenge"
meta_summary$group[meta_summary$Gene %in% Lysosome_pathways] <- "Lysosome pathways"
meta_summary$group[meta_summary$Gene %in% skin_melanin_enzymes] <- "Melanin enzymes"

# reorder genes
meta_summary_tmp <- meta_summary %>%
  group_by(Gene) %>%
  dplyr::summarise(MeanValue = mean(Avg))

meta_summary_tmp <- meta_summary_tmp[order(meta_summary_tmp$MeanValue,decreasing=T),]

# reorder the Genes within each group based on the mean
meta_summary$Gene <- factor(meta_summary$Gene, levels = meta_summary_tmp$Gene )

# get border colors of the sig unique DEG
meta_summary$Border <- "grey"
meta_summary$Border[paste0(meta_summary$Cell_Type,meta_summary$Gene) %in% paste0(res_uniq$Cell_type,res_uniq$Gene)] <- "red"

# get meta_summaries
meta_summary$Avg <- log10(meta_summary$Avg)

dplot <- meta_summary[meta_summary$group != "Pigmentation Network",]
dplot$Gene <- droplevels(dplot$Gene)
dp1 <- ggplot(dplot, aes(x=Gene, y=Cell_Type)) +
  geom_point(aes(size = Pct, fill = Avg, colour = Border), shape=21) +
  scale_size("% detected", range = c(0,6)) +
  scale_fill_gradientn(colours = viridisLite::mako(100),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Log10(Average\nexpression)") +
  scale_color_identity() + 
  ylab("") + xlab("") +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  facet_wrap(~group, scales = "free_x", nrow = 1)

dplot <- meta_summary[meta_summary$group == "Pigmentation Network",]
dplot$Gene <- droplevels(dplot$Gene)
dp2 <- ggplot(dplot, aes(x=Gene, y=Cell_Type)) +
  geom_point(aes(size = Pct, fill = Avg, colour = Border), shape=21) +
  scale_size("% detected", range = c(0,6)) +
  scale_fill_gradientn(colours = viridisLite::mako(100),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Log10(Average\nexpression)") +
  scale_color_identity() + 
  ylab("") + xlab("") +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  facet_wrap(~group, scales = "free_x", nrow = 1)


# save
# arrange <- ggarrange(plotlist=list(dp1,dp2), nrow=1, ncol=2)
ggsave("dotplot_kamath.ctr_yhgenes_1.pdf", dp1, width = 50, height = 20, units = "cm")
ggsave("dotplot_kamath.ctr_yhgenes_2.pdf", dp2, width = 50, height = 20, units = "cm")


#-----------  4)  Estimate the occurance of TYR expression by matching percentiles of UMI counts from TYR+/ TYR- ---------- 

# get counts
meta_dat$TYR <- exp_mat1["TYR",] > 0
group1_umi_counts <- meta_dat$nCount_RNA[meta_dat$Cell_Type == "SOX6_AGTR1" & meta_dat$TYR == "TRUE"]
group2_umi_counts <- meta_dat$nCount_RNA[meta_dat$Cell_Type == "SOX6_AGTR1" & meta_dat$TYR == "FALSE"]
# Sort the counts
sorted_group1 <- sort(group1_umi_counts)
sorted_group2 <- sort(group2_umi_counts)
# Calculate probabilities for each sample in Group 2
probabilities_group2 <- diff(c(0, quantiles_group1))
# Normalize probabilities to sum to 1
probabilities_group2 <- probabilities_group2 / sum(probabilities_group2)
# Remove the fewest samples from Group 2 to match distribution
subsetted_group2 <- sorted_group2[-sample(seq_along(sorted_group2), size = length(sorted_group2) - length(sorted_group1), replace = TRUE, prob = probabilities_group2)]
# predicted % of SOX6_AGTR1 cells with TYR expression.
length(group1_umi_counts)/(length(subsetted_group2) + length(group1_umi_counts)) * 100



#----------- 5) DEG analysis between TYR+/TYR- SOX6_AGTR1 cells of similar sequencing depth ----------- 

meta_dat$TYR <- exp_mat1["TYR",] > 0
tmp <- meta_dat[meta_dat$Cell_Type == "SOX6_AGTR1" & meta_dat$TYR == "TRUE",]
tmp2 <- meta_dat[meta_dat$Cell_Type == "SOX6_AGTR1" & meta_dat$TYR == "FALSE" & meta_dat$nCount_RNA > quantile(tmp$nCount_RNA)[4],]
t.test(tmp$nCount_RNA,tmp2$nCount_RNA)
meta_dat2 <- rbind(tmp,tmp2)
exp_mat2 <- as.matrix(seurat_obj[['RNA']]$counts[,row.names(meta_dat2)])


# format covariates
meta_dat2$TYR <- as.factor(meta_dat2$TYR)
meta_dat2$Cell_Type <- as.factor(meta_dat2$Cell_Type)
meta_dat2$donor_id <- as.factor(meta_dat2$donor_id)
meta_dat2$Donor_PMI<- as.factor(meta_dat2$Donor_PMI)
meta_dat2$libname <- as.factor(meta_dat2$libname)
meta_dat2$sex <- as.factor(meta_dat2$sex)
meta_dat2$nCount_RNA <- as.numeric(meta_dat2$nCount_RNA)
meta_dat2$Donor_Age <- as.numeric(meta_dat2$Donor_Age)

## Run Voom
design <- model.matrix( ~ TYR + donor_id + libname + Donor_Age + sex + nCount_RNA, data=meta_dat2,droplevels = F)
colnames(design) <- make.names(colnames(design))

# fit design
v <- voom(exp_mat2,design)
vfit <- lmFit(v)

# get results and assign to list
cont.matrix <- makeContrasts(A="TYRTRUE",levels=design)
fit2 <- contrasts.fit(vfit, cont.matrix)
vfit2 <- eBayes(fit2)
options(digits=3)

# Select significant upregulated DEGs and assign to list
tmp <- topTable(vfit2,number = Inf,sort.by="P")
res_tyr.sox6agtr1 <- tmp[tmp$adj.P.Val < 0.05 & abs(tmp$logFC) > 1,]

save(res_tyr.sox6agtr1,file="limma_res_tyr.sox6agtr1.Rdata")

## volcano plot highlighting specific genes
# keyvals <- ifelse(
#   row.names(lm_res) %in% dopamine, 'purple',
#   ifelse(row.names(lm_res) %in% ca_transmission, 'royalblue',
#          'grey'))
# keyvals[is.na(keyvals)] <- 'black'
# names(keyvals)[keyvals == 'purple'] <- 'dopamine'
# names(keyvals)[keyvals == 'royalblue'] <- 'Ca2+ & transmission'

lm_res <- tmp[tmp$adj.P.Val < 0.05,]
g1 <- EnhancedVolcano(lm_res,
                      lab = rownames(lm_res),
                      x = 'logFC',
                      y = 'adj.P.Val',
                      #selectLab = c(dopamine,ca_transmission),
                      xlab = bquote(~Log[2]~ 'fold change'),
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      #pointSize = 5,
                     # labSize = 4.5,
                      # shape = c(6, 4, 2, 11),
                      #colCustom = keyvals,
                      colAlpha = 1,
                      legendPosition = 'left',
                      legendLabSize = 15,
                      legendIconSize = 5.0,
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'black',
                      arrowheads = FALSE,
                      gridlines.major = TRUE,
                      gridlines.minor = FALSE,
                      border = 'partial',
                      borderWidth = 1.5,
                      borderColour = 'black') 
# save
pdf(paste0("volcano_sox6.agtr1.tyrposVtyrneg.pdf"))
g1
dev.off()


g2 <- EnhancedVolcano(lm_res,
                      lab = rownames(lm_res),
                      x = 'logFC',
                      y = 'adj.P.Val',
                      #selectLab = c(dopamine,ca_transmission),
                      xlab = bquote(~Log[2]~ 'fold change'),
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      #pointSize = 5,
                      # labSize = 4.5,
                      # shape = c(6, 4, 2, 11),
                      #colCustom = keyvals,
                      colAlpha = 1,
                      legendPosition = 'left',
                      legendLabSize = 15,
                      legendIconSize = 5.0,
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'black',
                      arrowheads = FALSE,
                      gridlines.major = TRUE,
                      gridlines.minor = FALSE,
                      border = 'partial',
                      borderWidth = 1.5,
                      borderColour = 'black') +ylim(0,40)

# save
pdf(paste0("volcano_sox6.agtr1.tyrposVtyrneg_zoom.pdf"))
g2
dev.off()


#----------- 6) violin plot of DEG interest -----------

# load subsetted and normalised control Kamath data
public_dataset_location = "/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/" # this is the mapped project drive
load(paste0(public_dataset_location,"kamath_seurat_sub.Rdata"))
exp_dat_sub <- as.matrix(seurat_obj_kamath_sub[['SCT']]$data)[row.names(seurat_obj_kamath_sub[['SCT']]$data) %in% genes_to_plot,]
meta_dat_sub <- as.data.frame(seurat_obj_kamath_sub@meta.data)
data_table <- cbind(meta_dat_sub,t(exp_dat_sub))
data_table <- data_table[,!duplicated(colnames(data_table))]
data_table <- data_table[ grep("SOX6_|CALB1_",meta_dat_sub$Cell_Type),]

gene_name <- c("MITF")
#data_table$log10_gene <- data_table[,gene_name]
x_variable = "Cell_Type"
y_variable = gene_name
x_lab = "Cell Type"
y_lab = "Norm. Counts"
colour_palette = rev(magma(30))[1:length(unique(dplot[,x_variable]))]

# format for plotting
tmp <- aggregate(formula(paste(gene_name," ~ ",x_variable)), data_table, mean)
fact_lvls <- tmp[order(-tmp[,gene_name]),][,x_variable]
data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)


bxp <- ggviolin(
  data_table, x = x_variable, y = y_variable, 
  fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 10)) + ylab(y_lab) + xlab(x_lab) + 
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black") + ylim(0,4)

# save plot
arrange <- ggarrange(plotlist=list(bxp), nrow=2, ncol=2, widths = c(2,2))
ggsave("Vln_kamath_MITF.CTR.DA.pdf", arrange, width = 8, height = 6)



#----------- 7) Run Gene Ontology -----------

exp_cont <- res_tyr.sox6agtr1
#exp_cont <- tmp[tmp$adj.P.Val < 0.05 & abs(tmp$logFC) > 0.6,]
exp_cont$Gene <- row.names(exp_cont)

# get unique
tmp <- exp_cont

# upreg list
upreg <- exp_cont$Gene[exp_cont$logFC > 0]
length(upreg)

# select databases
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015","Jensen_COMPARTMENTS")

# upreg top 3 passing adj p-val < 0.05
enriched <- enrichr(exp_cont$Gene, dbs)
lapply(enriched, function(x) {
  #x <- x[x$Adjusted.P.value < 0.05,]
  x[1:3,]
})

# plot enriched pathways
p1 <- plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value", title = "GO Molecular Function")
p2 <- plotEnrich(enriched[["GO_Cellular_Component_2015"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value", title = "GO Cellular Component")
p3 <- plotEnrich(enriched[["Jensen_COMPARTMENTS"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "Adjusted.P.value", title = "Jensen COMPARTMENTS")

arrange <- ggarrange(plotlist=list(p1,p2,p3), nrow=2, ncol=2, widths = c(2,2))
# ggsave("GO_enrichment_upreg_DAneuron.png", arrange)

ggsave("GO_Molecular_Function_2015_SOX6.AGTR1.TYRpos.png", p1)
ggsave("GO_Cellular_Component_2015_SOX6.AGTR1.TYRpos.png", p2)
ggsave("Jensen_COMPARTMENTS_SOX6.AGTR1.TYRpos.png", p3)

# write to excel
# write.xlsx(enriched[["GO_Molecular_Function_2015"]], file = "DEG_GO_Molecular_Function_2015_upreg_DAneuron6.12.xlsx", row.names = FALSE)
# write.xlsx(enriched[["GO_Cellular_Component_2015"]], file = "DEG_GO_Cellular_Component_2015_upreg_DAneuron6.12.xlsx", row.names = FALSE)
# write.xlsx(enriched[["Jensen_COMPARTMENTS"]], file = "DEG_Jensen_COMPARTMENTS_upreg_DAneuron6.12.xlsx", row.names = FALSE)

## test for genes interest
res <- tmp[tmp$adj.P.Val < 0.05,]
res[which(row.names(res_tyr.sox6agtr1) %in% yf_pigment_genes),]
res_tyr.sox6agtr1[which(row.names(res_tyr.sox6agtr1) %in% yf_pigment_genes),]

res_tyr.sox6agtr1[c("C1GALT1C1","GALNT5","GALNT13","CHPF2","B4GALNT1"),]




# ## heatmap
library(cluster)
# load(rdata)
# dat <- as.matrix(gxdat_s@assays$RNA$data)
# meta_dat$Diagnosis<- factor(meta_dat$Diagnosis, levels=c('CTR','ILBD','ePD','lPD'))

# subset for genes of interest
tmp3 <- tmp2[order(tmp2$adj.P.Val),]
genes_interest <- c(tmp3$Gene[1:10])
tmp <- dat[genes_interest,]

# subset for AOIs of interest
meta_small <- meta_dat[meta_dat$segment %in% "Full ROI" & meta_dat$ROI == "SND",]
dplot <- log10(tmp[,row.names(meta_small)])

# column annotations
anno_df = data.frame(
  Dx = meta_small$Diagnosis,
  Region = meta_small$ROI
)

ha = HeatmapAnnotation(df = anno_df)

# draw heatmap
#pdf("heatmap_snrnaseq_channelmembers_median_annot.pdf", width=15, height=15)
hm <- Heatmap(dplot, 
              top_annotation = ha,
              cluster_rows = diana(dplot)
              
              # cluster_columns = diana(t(dplot))
              #column_split=as.factor(Cell_Broad)
)
hm

draw(hm,
     column_title = "Genes interest in Kamath et al",
     column_title_gp=grid::gpar(fontsize=16))
#dev.off()













# 4) distribution of ncounts per gene highlighting TYR in LBD and PD
# get expression data and cell metadata
meta_dat <- meta_dat_total[meta_dat_total$Cell_Type %in% unique(meta_dat_total$Cell_Type)[grep("SOX6_|CALB1_",unique(meta_dat_total$Cell_Type))],]
exp_mat1 <- exp_dat[,row.names(meta_dat)]

# Create a contingency table and Perform Fisher's Exact Test
# Create a logical vector for the specific cell type
Cell_type <- unique(meta_dat$Cell_Type)
for (val in 1:length(Cell_type)){
  tmp <- meta_dat[meta_dat$Cell_Type == Cell_type[val],]
  tmp2 <- exp_mat1[,row.names(tmp)]
  x <- tmp$Status == "Ctrl"
  y <- tmp2["TYR", ] >= 1
  contingencyTable <- table(x,y)
  testResult <- fisher.test(contingencyTable)
  print(Cell_type[val])
  print(testResult)
}

# t-test on depth
x = meta_dat$nCount_RNA[meta_dat$Cell_Type == "SOX6_AGTR1" & meta_dat$Status == "Ctrl"]
y = meta_dat$nCount_RNA[meta_dat$Cell_Type == "SOX6_AGTR1" & meta_dat$Status %in% c("LBD","PD")]
t.test(x,y)

# Create boxplot TYR expression by cell-type and Dx
# Calculate the percentage of each cell type with counts > 0 for TYR gene
gene_expression_df <- as.data.frame(exp_mat1)
cell_types_vector <- meta_dat$Cell_Type
diagnosis_vector <- meta_dat$Status

tmp <- gene_expression_df %>%
  filter(row.names(exp_mat1) == "TYR") %>%
  pivot_longer(cols = everything(), names_to = "Cell", values_to = "Count") %>%
  mutate(HasCount = Count > 0) %>%
  group_by(CellType = cell_types_vector, Diagnosis = diagnosis_vector )

dplot <- tmp %>%
  group_by(CellType, Diagnosis) %>%
  dplyr::summarize(Percentage = mean(HasCount) * 100)

dplot$CellType <- factor(dplot$CellType, levels = unique(dplot$CellType[order(-dplot$Percentage)]))
custom_colors <- c("Ctrl" = "grey", "LBD" = "#00AFBB", "PD" = "red")

g3 <- ggplot(dplot,aes( x = CellType, y = Percentage, fill=Diagnosis)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percentage DA neurons > 0 TYR counts",
       x = "Cell Type",
       y = "% Cells > 0 TYR counts",
       fill = "Count > 0") +
  theme_minimal() + scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_fill_manual(values = custom_colors) + 
  geom_text(label="***", x=3,y=2) 


arrange <- ggarrange(plotlist=list(g3), nrow=2, ncol=2)
ggsave("plots_kamath.ctrpd.counts.pdf",arrange)





###########################################################
### Evaluate DA & NE neurons from Kamath & Webber et al ###
###########################################################

###########################################################
### 1) Format snRNAseq data from Kamath & Webber et al ###
###########################################################
# # 1) Kamath Control data
# public_dataset_location = "/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/" # this is the mapped project drive
# working_directory = "/Users/zacc/USyd/spatial_transcriptomics/analysis/public_datasets"
# setwd(working_directory)
# cell_ranger_data = paste0(public_dataset_location,"kamath_2023/") # single-cell cell ranger output
# # i) load single-cell data
# sc_obj <- Read10X(cell_ranger_data) # single-cell cell ranger output
# # ii) create seurat object
# seurat_obj <- CreateSeuratObject(sc_obj)
# rm(sc_obj)
# # iii) update metadata
# tsv <- read.delim(paste0(cell_ranger_data,"METADATA_PD.tsv")) # metadata
# tsv <- tsv[-1,]
# row.names(tsv) <- tsv$NAME
# seurat_obj <- AddMetaData(object = seurat_obj, metadata = tsv)
# umap_file <- list.files(cell_ranger_data, pattern="*UMAP.tsv", full.names=TRUE) # read in umap of cell IDs
# cell_id <- do.call(rbind, lapply(umap_file,fread))
# cell_id <- cell_id[!cell_id$NAME == "TYPE",]
# row.names(cell_id) <- cell_id$NAME
# seurat_obj <- AddMetaData(object = seurat_obj, metadata = cell_id)
# # iv) select control subjects
# seurat_obj_sub <- subset(seurat_obj, subset = Status == 'Ctrl'  )
# seurat_obj_sub <- subset(seurat_obj_sub, subset = organ__ontology_label == "substantia nigra pars compacta")
# # v) get DA neurons and normalize with SCT
# cell_types <- names(table(seurat_obj_sub@meta.data$Cell_Type))
# cell_types <- cell_types[grep("SOX6_|CALB1_",cell_types)]
# seurat_obj_sub <- subset(seurat_obj_sub, subset = Cell_Type %in% cell_types )
# #rm(seurat_obj)
# seurat_obj_sub <- SCTransform(seurat_obj_sub)
# save(seurat_obj_sub, file = "kamath_seurat_sub_DA.CTR.Rdata")
# 
# 
# # 2) Combine with Webber Control data
# load("/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/webber_re_seurat.Rdata")
# # i) merge datasets
# seurat_list <- c(seurat_obj_sub,seurat_webber)
# # ii) get variable features
# seurat_list <- lapply(X = seurat_list, function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)) 
# # iii) select features that are repeatedly variable across datasets for integration 
# features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000)
# seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
# 
# # iv) find anchors and intergrate data
# merge.anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = 'SCT', anchor.features = features)
# merge.combined.sct <- IntegrateData(anchorset = merge.anchors, normalization.method = 'SCT')
# 
# # v) Run the standard workflow for visualization and clustering
# merge.combined.sct <- ScaleData(merge.combined.sct, verbose = FALSE)
# merge.combined.sct <- RunPCA(merge.combined.sct, npcs = 30, verbose = FALSE)
# merge.combined.sct <- RunUMAP(merge.combined.sct, reduction = "pca", dims = 1:30)
# 
# merge.combined.sct@meta.data$dataset_merge <- c(rep("Kamath",15684),rep("Webber",20191))
# merge.combined.sct <- RunHarmony(merge.combined.sct, group.by.vars = c("dataset_merge") )
# merge.combined.sct <- FindNeighbors(merge.combined.sct, reduction = "harmony", dims = 1:30)
# merge.combined.sct <- FindClusters(merge.combined.sct, resolution = 0.5)
# 
# # vi) add metadata and save 
# merge.combined.sct@meta.data$cell_type_merge <- c(merge.combined.sct@meta.data$Cell_Type[!is.na(merge.combined.sct@meta.data$Cell_Type)],
#                                                   merge.combined.sct@meta.data$label_merged[is.na(merge.combined.sct@meta.data$Cell_Type)])
# 
# save(merge.combined.sct, file = "merge_kamath.webber_seurat_neurons.CTR.Rdata")

# # 3) Load data/ format for downstream analysis
# i) get raw count data and metadata for genes and cells of interest for linear modeling
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/public_datasets/merge_kamath.webber_seurat_neurons.CTR.Rdata")

# get cell-types for analysis
cell_types <- c( "CALB1_CALCR","CALB1_CRYM_CCDC68","CALB1_GEM","CALB1_PPP1R17",
                 "CALB1_RBP4","CALB1_TRHR","SOX6_AGTR1","SOX6_DDT",
                 "SOX6_GFRA2","SOX6_PART1","NE")

# ii) define count data for whole transcriptome DEG analysis
exp_dat1 <- as.data.frame(as.matrix(merge.combined.sct@assays$RNA$counts.SeuratProject))
exp_dat2 <- as.data.frame(as.matrix(merge.combined.sct@assays$RNA$counts.Siletti_recode))

common_genes <- intersect(row.names(exp_dat1),row.names(exp_dat2)) # ensure common genes
count_mat <- cbind(exp_dat1[common_genes,],exp_dat2[common_genes,]) # count matrix
plot_mat <- as.matrix(merge.combined.sct@assays$SCT$data[row.names(merge.combined.sct@assays$SCT$data) %in% common_genes,]) # plot matrix
common_genes <- intersect(row.names(count_mat),row.names(plot_mat)) 
count_mat <- count_mat[common_genes,]
plot_mat <- plot_mat[common_genes,]

meta_dat <- as.data.frame(merge.combined.sct@meta.data) # meta data

# select cells
count_mat <- count_mat[,meta_dat$cell_type_merge %in% cell_types ]
plot_mat <- plot_mat[,meta_dat$cell_type_merge %in% cell_types ]
meta_dat <- meta_dat[meta_dat$cell_type_merge %in% cell_types, ]

# checks
table(colnames(count_mat) == row.names(meta_dat))
table(row.names(count_mat) == row.names(plot_mat))
table(colnames(count_mat) == colnames(plot_mat))


#########################################################
### 2) Evaluate sequencing depth and TYR expression ###
#########################################################
# 1) Evaluate sequencing depth
aggregate(meta_dat$nCount_RNA,list(meta_dat$Cell_Type == "NE"), mean)

bxp <- ggboxplot(meta_dat, x="cell_type_merge", y="nCount_RNA",fill="cell_type_merge",palette = Cell_col,
  add = "mean_sd", size = 0.5) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10),
        plot.title = element_text(size = 10)) +
  ylab("Counts") + xlab("") +
  ggtitle("") 

arrange <- ggarrange(plotlist=list(bxp), nrow=2, ncol=2, widths = c(2,2))
ggsave("ncount_snRNAseq_webber.kamath_da.na.pdf", arrange,width = 8, height = 6)


# 2) Calculate the percentage of each cell type with counts > 0 for TYR gene
cell_types_vector <- meta_dat$cell_type_merge 
tmp <- count_mat %>%
  filter(row.names(count_mat) == "TYR") %>%
  pivot_longer(cols = everything(), names_to = "Cell", values_to = "Count") %>%
  mutate(HasCount = Count > 0) %>%
  group_by(CellType = cell_types_vector)

tmp2 <- table(tmp$HasCount,tmp$CellType)
dplot <- as.data.frame(cbind(Percentage = as.numeric(tmp2["TRUE",]/colSums(tmp2) * 100), CellType = colnames(tmp2)))
dplot$Percentage <- as.numeric(dplot$Percentage)
dplot$CellType <- factor(dplot$CellType, levels = dplot$CellType[order(-dplot$Percentage)])

g2 <- ggbarplot(dplot, x="CellType", y="Percentage",fill="CellType",palette = Cell_col,
        size = 0.5) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10),
        plot.title = element_text(size = 10)) +
  ylab("% Cells > 0 TYR counts") + xlab("") +
  ggtitle("")  + 
  geom_text(label="***", x=1,y=1.8) + ylim(0,2)

arrange <- ggarrange(plotlist=list(g2), nrow=2, ncol=2)
ggsave("plots_kamath.webber.ctr.counts.pdf",arrange,width = 8, height = 6)


#########################################################
### 3) Plot heatmaps of genes of interest ###
#########################################################
# select genes
genes_to_plot <- row.names(plot_mat)[row.names(plot_mat) %in% yf_pigment_genes] # total 125/ 130 genes interest

# # genes
# res_plot <- lapply(res_group, function(x) x[x$adj.P.Val < 0.05 & abs(x$logFC) > 1, ])
# genes_to_plot <- do.call(rbind,res_plot)
# genes_to_plot <- unique(genes_to_plot$Gene)

# define cells
y_list <- unique(meta_dat$cell_type_merge)

# aggregate data
x = t(plot_mat[genes_to_plot,])
y = meta_dat
dplot <- as.data.frame(x) %>%
  group_by(y$cell_type_merge) %>% 
  summarise_all("mean")

# # plot data frames
dplot1 <- t(dplot[,2:ncol(dplot)])
dplot1 <- scale(dplot1) # z-score
#dplot1 <- log10(dplot1) # log10

# remove genes not expressed (below zero z-score in all)
dplot1 <- dplot1[!apply(apply(dplot1,2,function(x) x < 0 ),1,all),]
  
# column annotations
anno_df = data.frame(
  Cell = dplot$`y$cell_type_merge`,
  Region = c("SND","SNV","VTA","SND",
             "VTA","VTA","LC","SNV","SNV",
             "VTA","SNV")
)

ha = HeatmapAnnotation(df = anno_df,
                       col = list(Cell = Cell_col,
                                  Region = ROI_col),
                       border = TRUE)

# row annotations
index <- which(row.names(dplot1) %in% genes_to_plot) 
annot_lab <- row.names(dplot1)[index]
row_ha = rowAnnotation(foo = anno_mark(at = index, 
                                   labels = annot_lab))

# row labels
tmp <- row.names(dplot1)
tmp[tmp %in%  other] <- "Genes interest"
tmp[tmp %in%  pigmentation_network] <- "Pigmentation Network"
tmp[tmp %in%  upstream_euNM] <- "Upstream euNM"
tmp[tmp %in%  c(CA_precursor_NM_stock_trans_metab_A9.A10,CA_precursor_NM_stock_trans_metab_A6)] <- "CA precursor"
tmp[tmp %in%  c(CA_functional_DA,CA_functional_NE)] <- "CA Functional"
tmp[tmp %in%  Stress_granules.free_radical_scavenging] <- "Stress Granule/ Scavenge"
tmp[tmp %in%  Lysosome_pathways] <- "Lysosome pathways"
tmp[tmp %in%  skin_melanin_enzymes] <- "Melanin enzymes"
tmp[tmp %in% yf_genes] <- "Genes interest"

# draw heatmap
pdf("heatmap_kamath.webber_n125.pdf",width = 10, height = 10)
hm <- Heatmap(dplot1, cluster_columns = TRUE,
              top_annotation = ha,
              row_split = tmp,
              col=viridis(100),
              #row_title_rot = 45,
              column_names_rot = 45,
              #right_annotation = row_ha,
              heatmap_legend_param = list(title = "z-score"),
              column_km = 4, row_names_gp = gpar(fontsize = 8),
              width = ncol(dplot1)*unit(5, "mm"), 
              height = nrow(dplot1)*unit(3, "mm")
)
draw(hm)
dev.off()


#############################################################################
### 4) DEG analysis between cell-types/ cell clusters / SOX6_AGTR1 TYR+/- ###
#############################################################################
# NOTES: 
# Many methods have been developed for DEG analysis in bulk RNA seq and single-cell RNAseq
# In a comprehensive comparison between methods it was shown that the edgeRQLFDetRate method performed best
# using a variety of performance metrics  - https://www.nature.com/articles/nmeth.4612
# however we compared this to limma voom and foud that expected DEGs within NE (DBH) were not positive
# we therefore used limma voom which also exhibit a good all round performance.

## i) DEGs between cell-type
# select genes
genes_to_plot <- row.names(count_mat)[row.names(count_mat) %in% yf_pigment_genes] # total 125/ 130 genes interest

# # format covariates
meta_dat$Cell_Type <- as.character(meta_dat$cell_type_merge)
meta_dat$Cell_Type <- as.factor(meta_dat$Cell_Type)
meta_dat$orig.ident <- as.factor(meta_dat$orig.ident)
meta_dat$nCount_RNA <- as.numeric(meta_dat$nCount_RNA)
meta_dat$dataset_merge <- as.factor(meta_dat$dataset_merge)

# limma voom
unique_cell_type <- unique(meta_dat$Cell_Type)
res <- list()
for (val in 1:length(unique_cell_type)){  
  print(unique_cell_type[val])
  meta_dat$condt <- meta_dat$Cell_Type == unique_cell_type[val]
 
  dge <- DGEList(count_mat[genes_to_plot,], group = meta_dat$condt)
  dge <- calcNormFactors(dge)
  design <- model.matrix( ~ condt, data=meta_dat)
  vm <- voom(dge, design = design, plot = TRUE)
  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)
  tt <- topTable(fit, n = Inf, adjust.method = "BH")
  res[[val]] <- tt
}

# name results
names(res) <- unique_cell_type
res <- lapply(res, function(x) {
  x$Gene <- row.names(x)
  row.names(x) <- NULL
  return(x)
})

res_cell <- res


## ii) DEGs between cell-type clusters
# define cell-type clusters - established from unsupervised heirarchical clustering from heatmaps
k2 <- c("CALB1_RBP4","SOX6_DDT")
k3 <- c("SOX6_AGTR1","CALB1_CRYM_CCDC68","SOX6_PART1","CALB1_PPP1R17")
k1 <- c("CALB1_GEM","CALB1_TRHR","SOX6_GFRA2","CALB1_CALCR")

# format covariates
meta_dat$Cell_Type <- as.character(meta_dat$cell_type_merge)
meta_dat$group <- meta_dat$Cell_Type
meta_dat$group[meta_dat$Cell_Type %in% k1 ] <- c("k1")
meta_dat$group[meta_dat$Cell_Type %in% k2 ] <- c("k2")
meta_dat$group[meta_dat$Cell_Type %in% k3 ] <- c("k3")
meta_dat$group <- as.factor(meta_dat$group)

# limma voom
unique_cell_type <- unique(meta_dat$group)
res <- list()
for (val in 1:length(unique_cell_type)){  
  print(unique_cell_type[val])
  meta_dat$condt <- meta_dat$group == unique_cell_type[val]
  
  dge <- DGEList(count_mat[genes_to_plot,], group = meta_dat$condt)
  dge <- calcNormFactors(dge)
  design <- model.matrix( ~ condt, data=meta_dat)
  vm <- voom(dge, design = design, plot = TRUE)
  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)
  tt <- topTable(fit, n = Inf, adjust.method = "BH")
  res[[val]] <- tt
}

# name results
names(res) <- unique_cell_type
res <- lapply(res, function(x) {
  x$Gene <- row.names(x)
  row.names(x) <- NULL
  return(x)
})

res_group <- res


## iii) DEGs between SOX6_AGTR1 TYR+ v TYR-
# # get total dataset
# exp_dat1 <- as.data.frame(as.matrix(merge.combined.sct@assays$RNA$counts.SeuratProject))
# exp_dat2 <- as.data.frame(as.matrix(merge.combined.sct@assays$RNA$counts.Siletti_recode))
# genes_to_plot <- intersect(row.names(exp_dat1),row.names(exp_dat2))
# exp_dat2 <- cbind(exp_dat1[genes_to_plot,],exp_dat2[genes_to_plot,])
# meta_dat_total <- as.data.frame(merge.combined.sct@meta.data)
# table(colnames(exp_dat2) == row.names(meta_dat_total))
# 
# # get merged count data for plotting
# plot_mat2 <- as.matrix(merge.combined.sct@assays$SCT$data[row.names(merge.combined.sct@assays$SCT$data) %in% genes_to_plot,])
# genes_to_plot <- intersect(row.names(exp_dat2), row.names(plot_mat2))
# exp_mat2 <- exp_dat2[genes_to_plot,]
# plot_mat2 <- plot_mat2[genes_to_plot,]

# define TYR- cells with high counts
# note; when selecting TYR- cells to contrast based on just counts, they have a lower feature/count ratio than TYR+ 
# indicating we are selecting potentially over-amplified samples. Hence select cells with high feature/count ratio & high features.

meta_dat_sub <- meta_dat[which(meta_dat$Cell_Type == "SOX6_AGTR1"),]
count_mat_sub <- count_mat[,row.names(meta_dat_sub)]
meta_dat_sub$condt <- c(count_mat_sub["TYR",] > 0)

tmp <- meta_dat_sub[meta_dat_sub$condt == "TRUE",]
tmp2 <- meta_dat_sub[meta_dat_sub$condt == "FALSE" & 
                       (meta_dat_sub$nFeature_RNA/meta_dat_sub$nCount_RNA) > quantile((tmp$nFeature_RNA/tmp$nCount_RNA))[3] & 
                       meta_dat_sub$nFeature_RNA > quantile(tmp$nFeature_RNA)[3],]

t.test(tmp$nFeature_RNA,tmp2$nFeature_RNA)
t.test((tmp$nFeature_RNA/tmp$nCount_RNA),
       (tmp2$nFeature_RNA/tmp2$nCount_RNA))

meta_dat$cell_type_tyr <- meta_dat$cell_type_merge
meta_dat$cell_type_tyr[row.names(meta_dat) %in% row.names(tmp)] <- "SOX6_AGTR1_TYRpos"
meta_dat$cell_type_tyr[row.names(meta_dat) %in% row.names(tmp2)] <- "SOX6_AGTR1_TYRneg"

# evaluate UMI and depth 
dplot <- meta_dat[meta_dat$cell_type_tyr %in% c("SOX6_AGTR1_TYRpos","SOX6_AGTR1_TYRneg"),]
p <- dplot %>%
  ggplot( aes(x=nCount_RNA, fill=cell_type_tyr)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
 # theme_ipsum() +
  labs(fill="")

p <- dplot %>%
  ggplot( aes(x=nFeature_RNA, fill=cell_type_tyr)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  # theme_ipsum() +
  labs(fill="")

dplot$feat.count <- dplot$nFeature_RNA/dplot$nCount_RNA

p <- dplot %>%
  ggplot( aes(x=feat.count, fill=cell_type_tyr)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  # theme_ipsum() +
  labs(fill="")

# limma voom
meta_dat_sub <- meta_dat[which(meta_dat$cell_type_tyr %in% c("SOX6_AGTR1_TYRpos","SOX6_AGTR1_TYRneg")),]
meta_dat_sub$condt <- meta_dat_sub$cell_type_tyr == "SOX6_AGTR1_TYRpos"
count_mat_sub <- count_mat[,row.names(meta_dat_sub)]

dge <- DGEList(count_mat_sub, group = meta_dat_sub$condt)
dge <- calcNormFactors(dge)
design <- model.matrix( ~ condt + sex, data=meta_dat_sub) # note: 101/109 sox6_agtr1 tyr+ cells are female
vm <- voom(dge, design = design, plot = TRUE)

vfit <- lmFit(vm)
colnames(design) <- make.names(colnames(design))
cont.matrix <- makeContrasts(A="condtTRUE",levels=design)
fit2 <- contrasts.fit(vfit, cont.matrix)
fit <- eBayes(fit2)
tt <- topTable(fit, n = Inf, adjust.method = "BH")

res_sox6agtr1.tyr <- tt


# fit <- lmFit(vm, design = design)
# fit <- eBayes(fit)
# tt <- topTable(fit, n = Inf, adjust.method = "BH")
# res_sox6agtr1.tyr <- tt

# save deg results
save(res_cell,res_group,res_sox6agtr1.tyr, file="deg_kamath.webber.Rdata")

# DEG cutoff
res_cut <- res_sox6agtr1.tyr[res_sox6agtr1.tyr$adj.P.Val < 0.05 & abs(res_sox6agtr1.tyr$logFC) > 0.5, ]
table(row.names(res_cut) %in% yf_pigment_genes)
row.names(res_cut)[row.names(res_cut) %in% yf_pigment_genes]

# write to file
write.table(res_cut, file = "deg_sox6agtr1_tyrposneg.txt", sep="\t", quote = F)

################################
### 5) Volcano plots of DEGs ###
################################
# volcano plotting function
volcano_plot <- function(lm_res) {
  vp1 <- EnhancedVolcano(lm_res,
                         lab = rownames(lm_res),
                         x = 'logFC',
                         y = 'adj.P.Val',
                         selectLab = c(yf_pigment_genes[yf_pigment_genes %in% row.names(res_cut)]),
                         xlab = bquote(~Log[2]~ 'fold change'),
                         shape = 16,
                         pCutoff = 0.05,
                         FCcutoff = 0.5,
                         pointSize = 5,
                         labSize = 4.5,
                         #ylim = c(0,20),
                         xlim = c(-2,2),
                         # shape = c(6, 4, 2, 11),
                         #colCustom = keyvals,
                         col = c('black', 'pink', 'purple', 'red3'),
                         colAlpha = 0.5,
                         legendPosition = 'none',
                         legendLabSize = 10,
                         legendIconSize = 5.0,
                         drawConnectors = TRUE,
                         widthConnectors = 0.5,
                         colConnectors = "grey30",
                         arrowheads = FALSE,
                         gridlines.major = TRUE,
                         gridlines.minor = FALSE,
                         border = 'partial',
                         borderWidth = 1.5,
                         borderColour = 'black') 
  return(vp1)
}

# plot
vp1 <- volcano_plot(res_sox6agtr1.tyr)

# save
pdf(paste0("volcano_sox6.agtr1.TYRposVTYRneg.pdf"))
vp1
dev.off()

#########################################
### 6) Gene Ontology analysis of DEGs ###
#########################################
# define upreg list
upreg <- res_cut[res_cut$logFC > 0,]
upreg$Gene <- row.names(upreg)
nrow(upreg)

# select databases
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015","Jensen_COMPARTMENTS")

# upreg top 10 passing adj p-val < 0.05
enriched <- enrichr(upreg$Gene, dbs)
enriched <- lapply(enriched, function(x) {
  x <- x[x$Adjusted.P.value < 0.05,]
  #x[1:10,]
})

# plot enriched pathways
p1 <- plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 10, numChar = 40, y = "Count", orderBy = "Adjusted.P.value", title = "GO Molecular Function")
p2 <- plotEnrich(enriched[["GO_Cellular_Component_2015"]], showTerms = 10, numChar = 40, y = "Count", orderBy = "Adjusted.P.value", title = "GO Cellular Component")
#p3 <- plotEnrich(enriched[["GO_Biological_Process_2015"]], showTerms = 10, numChar = 40, y = "Count", orderBy = "Adjusted.P.value", title = "GO Biological Process")
#p4<- plotEnrich(enriched[["Jensen_COMPARTMENTS"]], showTerms = 10, numChar = 40, y = "Count", orderBy = "Adjusted.P.value", title = "Jensen COMPARTMENTS")

# plot
arrange <- ggarrange(plotlist=list(p2,p1), nrow=3, ncol=1, widths = c(2,2))
pdf("GO_enrichment_downreg_sox6agtr1_TYR.pdf",width = 8, height = 6)
arrange
dev.off()


########################################################################################
### 7) Evaluate the SOX6_AGTR1 TYR+ genes of interest within NA cells and SOX6_AGTR1 TYR+ and TYR- ###
########################################################################################
#genes_to_plot <- c("EDN3","TYRP1","COA7","GSTM1","ITGB3","C4A","C4B")

# 1) evaluate all DEGs by aggregate in SOX6_AGTR1 and NA cells
genes_to_plot <- c(row.names(res_cut)[res_cut$logFC > 0.5],row.names(res_cut)[res_cut$logFC < -0.5])

# aggregate data
y = meta_dat[which(meta_dat$cell_type_tyr %in% c("SOX6_AGTR1_TYRpos","SOX6_AGTR1_TYRneg","NE")),]
x = t(plot_mat[genes_to_plot,row.names(y)])
dplot <- as.data.frame(x) %>%
  group_by(y$cell_type_tyr) %>% 
  summarise_all("mean")

# # plot data frames
dplot1 <- t(dplot[,2:ncol(dplot)])
dplot1 <- scale(dplot1)

# test correlation
cor.test(dplot1[,which(dplot$`y$cell_type_tyr` == "NE")],
         dplot1[,which(dplot$`y$cell_type_tyr` == "SOX6_AGTR1_TYRpos")])
cor.test(dplot1[,which(dplot$`y$cell_type_tyr` == "NE")],
         dplot1[,which(dplot$`y$cell_type_tyr` == "SOX6_AGTR1_TYRneg")])


# column annotations
anno_df = data.frame(
  Cell = dplot$`y$cell_type_tyr`)


ha = HeatmapAnnotation(df = anno_df,
                       col = list(Cell = c(Cell_col,Cell_col2)),
                       border = TRUE)

# draw heatmap
pdf("heatmap_kamath.webber_sx6agtr1_na_tyrdeg.pdf",width = 5, height = 5)
hm <- Heatmap(dplot1, cluster_columns = TRUE,
              top_annotation = ha,
              #row_split = tmp,
              col=viridis(100)[50:100],
              #row_title_rot = 45,
              column_names_rot = 45,
              #right_annotation = row_ha,
              heatmap_legend_param = list(title = "z-score"),
              row_km = 3, row_names_gp = gpar(fontsize = 4),
              width = ncol(dplot1)*unit(8, "mm"), 
              height = nrow(dplot1)*unit(0.2, "mm")
)
draw(hm)
dev.off()

# 2) evaluate up-regulated DEGs in NA cells
genes_to_plot <- c(row.names(res_cut)[res_cut$logFC > 0.5])

# get k-means clusters
y = meta_dat[which(meta_dat$cell_type_tyr %in% c("NE")),]
x = t(plot_mat[genes_to_plot,row.names(y)])
dplot1 <- x
cluster <- kmeans(dplot1, 2)
       

cluster$cluster

# assign cluster name
meta_dat$cell_type_tyr[row.names(meta_dat) %in% names(cluster$cluster)[cluster$cluster == 1] ] <- "NE_cluster1"
meta_dat$cell_type_tyr[row.names(meta_dat) %in% names(cluster$cluster)[cluster$cluster == 2] ] <- "NE_cluster2"


# boxplot of gene expression
gene_name <- c("ALCAM")
x_variable = "Cell_Type"
y_variable = gene_name
x_lab = "Cell Type"
y_lab = "Norm. Counts"
colour_palette = Cell_col

# format for plotting
tmp <- aggregate(formula(paste(gene_name," ~ ",x_variable)), data_table, mean)
fact_lvls <- tmp[order(-tmp[,gene_name]),][,x_variable]
data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)





gene_name <- c("ALCAM")

y = meta_dat[which(meta_dat$cell_type_tyr %in% c("SOX6_AGTR1_TYRpos","SOX6_AGTR1_TYRneg","NE_cluster1","NE_cluster2")),]
gene = t(plot_mat[gene_name,row.names(y)])

data_table <- cbind(y,gene)
data_table$gene <- as.numeric(data_table$gene)


bxp <- ggviolin(
  data_table, x = "cell_type_tyr", y = "gene", 
  fill = "cell_type_tyr"#, palette = colour_palette
  ) + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 10)) + 
  #ylab(y_lab) + xlab(x_lab) + 
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black") #+ ylim(1.09,1.12) + xlab("")

# save plot
arrange <- ggarrange(plotlist=list(bxp), nrow=2, ncol=2, widths = c(2,2))
ggsave("Vln_kamath.webber_MITF.CTR.DA.pdf", arrange, width = 8, height = 6)




x = t(plot_mat[genes_to_plot,row.names(y)])
dplot <- as.data.frame(x) %>%
  group_by(y$cell_type_tyr) %>% 
  summarise_all("mean")

# # plot data frames
dplot1 <- t(dplot[,2:ncol(dplot)])
dplot1 <- scale(dplot1)
#dplot1 <- log10(dplot1)

# column annotations
anno_df = data.frame(
  Cell = dplot$`y$cell_type_tyr`)


ha = HeatmapAnnotation(df = anno_df,
                       col = list(Cell = c(Cell_col,Cell_col2)),
                       border = TRUE)



dplot1 <- t(x)
# column annotations
anno_df = data.frame(
  Cell = y$cell_type_tyr)


ha = HeatmapAnnotation(df = anno_df,
                       col = list(Cell = c(Cell_col,Cell_col2)),
                       border = TRUE)




# row annotations
# index <- which(row.names(dplot1) %in% genes_to_plot) 
# annot_lab <- row.names(dplot1)[index]
# row_ha = rowAnnotation(foo = anno_mark(at = index, 
#                                        labels = annot_lab))

# row labels
# tmp <- row.names(dplot1)
# tmp[tmp %in%  other] <- "Genes interest"
# tmp[tmp %in%  pigmentation_network] <- "Pigmentation Network"
# tmp[tmp %in%  upstream_euNM] <- "Upstream euNM"
# tmp[tmp %in%  c(CA_precursor_NM_stock_trans_metab_A9.A10,CA_precursor_NM_stock_trans_metab_A6)] <- "CA precursor"
# tmp[tmp %in%  c(CA_functional_DA,CA_functional_NE)] <- "CA Functional"
# tmp[tmp %in%  Stress_granules.free_radical_scavenging] <- "Stress Granule/ Scavenge"
# tmp[tmp %in%  Lysosome_pathways] <- "Lysosome pathways"
# tmp[tmp %in%  skin_melanin_enzymes] <- "Melanin enzymes"
# tmp[tmp %in% yf_genes] <- "Genes interest"

# draw heatmap
pdf("heatmap_kamath.webber_sx6agtr1_tyrdeg.pdf",width = 100, height = 100)
hm <- Heatmap(dplot1, cluster_columns = TRUE,
              top_annotation = ha,
              #row_split = tmp,
              col=viridis(100)[50:100],
              #row_title_rot = 45,
              column_names_rot = 45,
              #right_annotation = row_ha,
              heatmap_legend_param = list(title = "z-score"),
              column_km = 2, row_names_gp = gpar(fontsize = 4),
              width = ncol(dplot1)*unit(1, "mm"), 
              height = nrow(dplot1)*unit(1, "mm")
)
draw(hm)
dev.off()














## 3) Select Genes of interest from 130 and plot heatmaps
# select genes above thresholds
res_plot <- lapply(res_group, function(x) x[x$adj.P.Val < 0.05 & abs(x$logFC) > 1, ])
genes_to_plot <- do.call(rbind,res_plot)
genes_to_plot <- unique(genes_to_plot$Gene)

library(ComplexHeatmap)

# cells
y_list <- unique(meta_dat$cell_type_merge)

# aggregate data
x = t(plot_mat)
y = meta_dat
dplot <- as.data.frame(x) %>%
  group_by(y$cell_type_merge) %>% 
  summarise_all("mean")

# # plot data frames
dplot1 <- t(dplot[,2:ncol(dplot)])
#dplot1 <- scale(dplot1)
#dplot1 <- log10(dplot1 + 0.5)
#dplot1 <- dplot1[genes_to_plot,]
colnames(dplot1) <- dplot$`y$cell_type_merge`

# column annotations
anno_df = data.frame(
  Cell = colnames(dplot1)
)

ha = HeatmapAnnotation(df = anno_df,
                       col = list(Cell = Cell_col))

# row annotations
index <- which(row.names(dplot1) %in% genes_to_plot) 
annot_lab <- row.names(dplot1)[index]
row_ha = rowAnnotation(foo = anno_mark(at = index, 
                                   labels = annot_lab))


# row labels
tmp <- row.names(dplot1)
tmp[tmp %in%  other] <- "Genes interest"
tmp[tmp %in%  pigmentation_network] <- "Pigmentation Network"
tmp[tmp %in%  upstream_euNM] <- "Upstream euNM"
tmp[tmp %in%  c(CA_precursor_NM_stock_trans_metab_A9.A10,CA_precursor_NM_stock_trans_metab_A6)] <- "CA precursor"
tmp[tmp %in%  c(CA_functional_DA,CA_functional_NE)] <- "CA Functional"
tmp[tmp %in%  Stress_granules.free_radical_scavenging] <- "Stress Granule/ Scavenge"
tmp[tmp %in%  Lysosome_pathways] <- "Lysosome pathways"
tmp[tmp %in%  skin_melanin_enzymes] <- "Melanin enzymes"
tmp[tmp %in% yf_genes] <- "Genes interest"

# draw heatmap
pdf("heatmap_kamath.webber_n130.n21.pdf",width = 10, height = 10)
hm <- Heatmap(dplot1, cluster_columns = TRUE,
              top_annotation = ha,
              row_split = tmp,
              #row_title_rot = 45,
              column_names_rot = 45,
              right_annotation = row_ha,
              heatmap_legend_param = list(title = "Mean Norm. Counts"),
              column_km = 4, row_names_gp = gpar(fontsize = 4),
              width = ncol(dplot1)*unit(5, "mm"), 
              height = nrow(dplot1)*unit(1.2, "mm")
)
draw(hm)
dev.off()









### violin plot of genes of interest
# load subsetted and normalised control Kamath data
public_dataset_location = "/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/" # this is the mapped project drive
load(paste0(public_dataset_location,"kamath_seurat_sub.Rdata"))
exp_dat_sub <- as.matrix(seurat_obj_kamath_sub[['SCT']]$data)[row.names(seurat_obj_kamath_sub[['SCT']]$data) %in% genes_to_plot,]
meta_dat_sub <- as.data.frame(seurat_obj_kamath_sub@meta.data)
data_table <- cbind(meta_dat_sub,t(exp_dat_sub))
data_table <- data_table[,!duplicated(colnames(data_table))]

data_table <- cbind(meta_dat,t(plot_mat))

gene_name <- c("MITF")
#data_table$log10_gene <- data_table[,gene_name]
x_variable = "Cell_Type"
y_variable = gene_name
x_lab = "Cell Type"
y_lab = "Norm. Counts"
colour_palette = Cell_col

# format for plotting
tmp <- aggregate(formula(paste(gene_name," ~ ",x_variable)), data_table, mean)
fact_lvls <- tmp[order(-tmp[,gene_name]),][,x_variable]
data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)


bxp <- ggboxplot(
  data_table, x = x_variable, y = y_variable, 
  fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 10)) + ylab(y_lab) + xlab(x_lab) + 
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="black") + ylim(0,5) + xlab("")

# save plot
arrange <- ggarrange(plotlist=list(bxp), nrow=2, ncol=2, widths = c(2,2))
ggsave("Vln_kamath.webber_MITF.CTR.DA.pdf", arrange, width = 8, height = 6)







# construct design matrix
design <- model.matrix( ~ 0 + Cell_Type + orig.ident + dataset_merge + nCount_RNA, data=meta_dat)
colnames(design) <- make.names(colnames(design))

# fit design
v <- voom(exp_mat,design)
vfit <- lmFit(v)

# run through non-baseline covariates
for (val in 1:length(unique(meta_dat$Cell_Type))){
  # Perform LIMMA contrasts
  cont.matrix <- makeContrasts(A=colnames(design)[val],levels=design)
  #cont.matrix <- makeContrasts(A="(groupk1 + groupk2 + groupk3) - groupNE",levels=design)
  fit2 <- contrasts.fit(vfit, cont.matrix)
  vfit2 <- eBayes(fit2)
  options(digits=3)
  
  # Select significant upregulated DEGs and assign to list
  tmp <- topTable(vfit2,number = Inf,sort.by="P")
  res[[val]] <- tmp
}

# name results
names(res) <- levels(meta_dat$group)
res <- lapply(res, function(x) {
  x$Gene <- row.names(x)
  row.names(x) <- NULL
  return(x)
})





# 4) perform DEG analysis between each cell-type
# # format covariates
meta_dat$Cell_Type <- as.character(meta_dat$cell_type_merge)
meta_dat$Cell_Type <- as.factor(meta_dat$Cell_Type)
meta_dat$orig.ident <- as.factor(meta_dat$orig.ident)
meta_dat$nCount_RNA <- as.numeric(meta_dat$nCount_RNA)
meta_dat$dataset_merge <- as.factor(meta_dat$dataset_merge)

## Run Voom
# 1) run first co variate by enabling accession in design matrix
meta_dat$Cell_Type <- factor(meta_dat$Cell_Type, levels = rev(levels(meta_dat$Cell_Type)))

design <- model.matrix( ~ Cell_Type + orig.ident + dataset_merge + nCount_RNA, data=meta_dat,droplevels = F)
colnames(design) <- make.names(colnames(design))

# fit design
v <- voom(exp_mat,design)
vfit <- lmFit(v)

# get results and assign to list
res <- list()
cont.matrix <- makeContrasts(A=colnames(design)[length(cell_types)],levels=design)
fit2 <- contrasts.fit(vfit, cont.matrix)
vfit2 <- eBayes(fit2)
options(digits=3)

# Select significant upregulated DEGs and assign to list
tmp <- topTable(vfit2,number = Inf,sort.by="P")
res[[1]] <- tmp

# 2) run rest of co variate by enabling accession in design matrix
meta_dat$Cell_Type <- factor(meta_dat$Cell_Type, levels = rev(levels(meta_dat$Cell_Type)))
design <- model.matrix( ~ Cell_Type + orig.ident + dataset_merge + nCount_RNA, data=meta_dat,droplevels = F)
colnames(design) <- make.names(colnames(design))

# fit design
v <- voom(exp_mat,design)
vfit <- lmFit(v)

# run through non-baseline covariates
for (val in 2:length(cell_types)){
  # Perform LIMMA contrasts
  cont.matrix <- makeContrasts(A=colnames(design)[val],levels=design)
  fit2 <- contrasts.fit(vfit, cont.matrix)
  vfit2 <- eBayes(fit2)
  options(digits=3)
  
  # Select significant upregulated DEGs and assign to list
  tmp <- topTable(vfit2,number = Inf,sort.by="P")
  res[[val]] <- tmp
}

# name results
names(res) <- levels(meta_dat$Cell_Type)
res <- lapply(res, function(x) {
  x$Gene <- row.names(x)
  row.names(x) <- NULL
  return(x)
})

res_cell <- res

# write to file
library(writexl)
write_xlsx(res_cell, path ="Kamath.Webber_ctr_DEG_pigment.n130.xlsx")




# 5) perform DEG analysis between each k-means group
res <- list()
# define k-means groups from heatmap
k2 <- c("CALB1_RBP4","SOX6_DDT")
k3 <- c("SOX6_AGTR1","CALB1_CRYM_CCDC68","SOX6_PART1","CALB1_PPP1R17")
k1 <- c("CALB1_GEM","CALB1_TRHR","SOX6_GFRA2","CALB1_CALCR")

# # format covariates
meta_dat$Cell_Type <- as.character(meta_dat$cell_type_merge)
meta_dat$group <- meta_dat$Cell_Type
meta_dat$group[meta_dat$Cell_Type %in% k1 ] <- c("k1")
meta_dat$group[meta_dat$Cell_Type %in% k2 ] <- c("k2")
meta_dat$group[meta_dat$Cell_Type %in% k3 ] <- c("k3")
meta_dat$group <- as.factor(meta_dat$group)
groups <- unique(meta_dat$group)

## Run Voom
# 1) run first co variate by enabling accession in design matrix
design <- model.matrix(~ 0 + group + orig.ident + dataset_merge + nCount_RNA, data=meta_dat)
colnames(design) <- make.names(colnames(design))

# fit design
v <- voom(exp_mat,design)
vfit <- lmFit(v)

# run through non-baseline covariates
for (val in 1:length(groups)){
  # Perform LIMMA contrasts
  cont.matrix <- makeContrasts(A=colnames(design)[val],levels=design)
  
  cont.matrix <- makeContrasts(A="(groupk1 + groupk2 + groupk3) - groupNE",levels=design)
  fit2 <- contrasts.fit(vfit, cont.matrix)
  vfit2 <- eBayes(fit2)
  options(digits=3)
  
  # Select significant upregulated DEGs and assign to list
  tmp <- topTable(vfit2,number = Inf,sort.by="P")
  res[[val]] <- tmp
}

# name results
names(res) <- levels(meta_dat$group)
res <- lapply(res, function(x) {
  x$Gene <- row.names(x)
  row.names(x) <- NULL
  return(x)
})






suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))


    dge <- DGEList(exp_mat, group = meta_dat$group)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~ 0 + group, data=meta_dat)
    y <- new("EList")
    y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
    fit <- lmFit(y, design = design)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    tt <- topTable(fit, n = Inf, adjust.method = "BH")
  
  hist(tt$P.Value, 50)
  hist(tt$adj.P.Val, 50)
  limma::plotMDS(dge, col = as.numeric(as.factor(meta_dat$group)), pch = 19)
  plotMD(fit)
  
  list(session_info = session_info,
       timing = timing,
       tt = tt,
       df = data.frame(pval = tt$P.Value,
                       padj = tt$adj.P.Val,
                       row.names = rownames(tt)))
}








## Stouffer
# Function to combine p-values using Stouffer's Method
combine_pvalues_stouffer <- function(p_values) {
  # Convert p-values to Z-scores
  z_scores <- qnorm(p_values / 2, lower.tail = FALSE)
  
  # Calculate combined Z-score
  combined_z <- sum(z_scores) / sqrt(length(z_scores))
  
  # Convert combined Z-score back to a p-value
  combined_p <- 2 * pnorm(-abs(combined_z))
  
  return(combined_p)
}

# get representative DEG for names
res_st <- res[[1]]

# calc stouffer
stouffer_p <- list()
for (z in 1:nrow(res_st)){
  gene <- res_st$Gene[z]
  p_values <- unlist(lapply(res,function(x) x[x[,7] == gene,][5]))
  stouffer_p[[z]] <- combine_pvalues_stouffer(p_values)
}

res_combined$stouffer_p <- combine_pvalues_stouffer(res_combined$adj.P.Val)









## 4) Plot heatmap of gene
library(ComplexHeatmap)

# cells
y_list <- unique(meta_dat$cell_type_merge)

# aggregate data
x = t(plot_mat)
y = meta_dat
dplot <- as.data.frame(x) %>%
  group_by(y$cell_type_merge) %>% 
  summarise_all("mean")

# # plot data frames
dplot1 <- t(dplot[,2:ncol(dplot)])
dplot1 <- scale(dplot1)
#dplot1 <- log10(dplot1)

# column annotations
anno_df = data.frame(
  Cell = dplot$`y$cell_type_merge`
)

ha = HeatmapAnnotation(df = anno_df,
                       col = list(Cell = Cell_col
                       ))
# row colors
tmp <- row.names(dplot1)
#tmp[tmp %in% yf_genes] <- "Genes interest"
tmp[tmp %in%  other] <- "Genes interest"
tmp[tmp %in%  pigmentation_network] <- "Pigmentation Network"
tmp[tmp %in%  upstream_euNM] <- "Upstream euNM"
tmp[tmp %in%  c(CA_precursor_NM_stock_trans_metab_A9.A10,CA_precursor_NM_stock_trans_metab_A6)] <- "CA precursor"
tmp[tmp %in%  c(CA_functional_DA,CA_functional_NE)] <- "CA Functional"
tmp[tmp %in%  Stress_granules.free_radical_scavenging] <- "Stress Granule/ Scavenge"
tmp[tmp %in%  Lysosome_pathways] <- "Lysosome pathways"
tmp[tmp %in%  skin_melanin_enzymes] <- "Melanin enzymes"
tmp[tmp %in% yf_genes] <- "Genes interest"


# m = matrix(rnorm(1000), nrow = 100)
# rownames(m) = 1:100
# ha = rowAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 97:100), 
#                                    labels = month.name[1:10]))
# Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha,
#         row_names_side = "left", row_names_gp = gpar(fontsize = 4))

# draw heatmap
pdf("heatmap_kamath.webber_n130.pdf",width = 10, height = 10)
hm <- Heatmap(dplot1, cluster_columns = TRUE,
              top_annotation = ha,
              row_split = tmp,
              heatmap_legend_param = list(title = "Log10(mean(Cell Prop.))"),
              column_km = 4, row_names_gp = gpar(fontsize = 4),
              width = ncol(dplot1)*unit(5, "mm"), 
              height = nrow(dplot1)*unit(1.2, "mm")
              )
draw(hm)
dev.off()







# # define any DEGs that are unique to a cell-type
# check <- 1:length(res)
# res_uniq <- list()
# for (i in 1:length(check)){
#   check_min <- check[check != i]
#   tmp <- do.call(rbind,res[check_min])
#   tmp2 <- as.data.frame(res[[i]])
#   
#   # unique to list
#   tmpa <- tmp[tmp$logFC > 0,]
#   tmp2a <- tmp2[tmp2$logFC > 0,]
#   tmp3 <- tmp2a[!tmp2a$Gene %in% tmpa$Gene,]
# 
#   tmpa <- tmp[tmp$logFC < 0,]
#   tmp2a <- tmp2[tmp2$logFC < 0,]
#   tmp4 <- tmp2a[!tmp2a$Gene %in% tmpa$Gene,]
# 
#   
#   res_uniq[[i]] <- rbind(tmp3,tmp4)
# }
# names(res_uniq) <- levels(meta_dat$Cell_Type)
# res_uniq <- do.call(rbind,res_uniq)
# res_uniq$Cell_type <-  unlist(lapply(strsplit(row.names(res_uniq),"\\."), function(x) x[1]))

## OR define any DEGs that are sig to cell-type
res_uniq <- do.call(rbind,res)
res_uniq <- res_uniq[res_uniq$adj.P.Val < quantile(res_uniq$adj.P.Val)[3] & res_uniq$logFC > quantile(res_uniq$logFC[res_uniq$logFC > 0])[3] ,]
res_uniq$Cell_type <-  unlist(lapply(strsplit(row.names(res_uniq),"\\."), function(x) x[1]))

# calculate percentile of each gene
total_expression <- rowSums(plot_mat)
percentile_exp <- ecdf(total_expression)(total_expression) * 100
names(percentile_exp) <- row.names(plot_mat)
sort(percentile_exp[names(percentile_exp) %in% genes_to_plot ])

# calculate the number of cells detected in
num_cells_expressed <- apply(plot_mat, 1, function(x) sum(x > 0))
names(num_cells_expressed) <- row.names(plot_mat)
sort(num_cells_expressed[names(num_cells_expressed) %in% genes_to_plot ])

# merge
meta1 <- meta_dat
meta1 <- bind_cols(meta1, as.data.frame(t(plot_mat)))

# select cols interest
meta1 <- meta1[,c("Cell_Type",genes_to_plot)]

# long format
meta <- pivot_longer(meta1, -Cell_Type, names_to="Gene", values_to="Expression")

# create summary
meta_summary <- meta %>%
  dplyr::group_by(Cell_Type, Gene) %>%
  dplyr::summarise(Avg = mean(Expression),
                   Pct = sum(Expression > 0) / length(Expression) * 100)

# create gene groups
meta_summary$group <- as.character(meta_summary$Gene)
meta_summary$group[meta_summary$Gene %in% yf_genes] <- "Genes interest"
meta_summary$group[meta_summary$Gene %in% other] <- "Genes interest"
meta_summary$group[meta_summary$Gene %in% pigmentation_network] <- "Pigmentation Network"
meta_summary$group[meta_summary$Gene %in% upstream_euNM] <- "Upstream euNM"
meta_summary$group[meta_summary$Gene %in% c(CA_precursor_NM_stock_trans_metab_A9.A10,CA_precursor_NM_stock_trans_metab_A6)] <- "CA precursor"
meta_summary$group[meta_summary$Gene %in% c(CA_functional_DA,CA_functional_NE)] <- "CA Functional"
meta_summary$group[meta_summary$Gene %in% Stress_granules.free_radical_scavenging] <- "Stress Granule/ Scavenge"
meta_summary$group[meta_summary$Gene %in% Lysosome_pathways] <- "Lysosome pathways"
meta_summary$group[meta_summary$Gene %in% skin_melanin_enzymes] <- "Melanin enzymes"

# reorder genes
meta_summary_tmp <- meta_summary %>%
  group_by(Gene) %>%
  dplyr::summarise(MeanValue = mean(Avg))

meta_summary_tmp <- meta_summary_tmp[order(meta_summary_tmp$MeanValue,decreasing=T),]

# reorder the Genes within each group based on the mean
meta_summary$Gene <- factor(meta_summary$Gene, levels = meta_summary_tmp$Gene )

# get border colors of the sig unique DEG
meta_summary$Border <- "grey"
meta_summary$Border[paste0(meta_summary$Cell_Type,meta_summary$Gene) %in% paste0(res_uniq$Cell_type,res_uniq$Gene)] <- "red"

# get meta_summaries
meta_summary$Avg <- log10(meta_summary$Avg)

dplot <- meta_summary[meta_summary$group != "Pigmentation Network",]
dplot$Gene <- droplevels(dplot$Gene)
dp1 <- ggplot(dplot, aes(x=Gene, y=Cell_Type)) +
  geom_point(aes(size = Pct, fill = Avg, colour = Border), shape=21) +
  scale_size("% detected", range = c(0,6)) +
  scale_fill_gradientn(colours = viridisLite::mako(100),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Log10(Average\nexpression)") +
  scale_color_identity() + 
  ylab("") + xlab("") +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  facet_wrap(~group, scales = "free_x", nrow = 1)

dplot <- meta_summary[meta_summary$group == "Pigmentation Network",]
dplot$Gene <- droplevels(dplot$Gene)
dp2 <- ggplot(dplot, aes(x=Gene, y=Cell_Type)) +
  geom_point(aes(size = Pct, fill = Avg, colour = Border), shape=21) +
  scale_size("% detected", range = c(0,6)) +
  scale_fill_gradientn(colours = viridisLite::mako(100),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Log10(Average\nexpression)") +
  scale_color_identity() + 
  ylab("") + xlab("") +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  facet_wrap(~group, scales = "free_x", nrow = 1)


# save
# # arrange <- ggarrange(plotlist=list(dp1,dp2), nrow=1, ncol=2)
ggsave("dotplot_kamath.webber.ctr_sctdata_yhgenes_1.pdf", dp1, width = 50, height = 20, units = "cm")
ggsave("dotplot_kamath.webber.ctr_sctdata_yhgenes_2.pdf", dp2, width = 50, height = 20, units = "cm")













