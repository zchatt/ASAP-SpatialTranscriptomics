# Details: Scripts used for human validation of DEGs identified by CHCHD2 point mutant (PM) mouse model, replicating this autosomal dominant form of PD. 

library(readxl)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(Seurat)
library(SeuratObject)
library(viridis)
library(dplyr)
library(ggpubr)
library(stargazer)
library(edgeR)

############################################################################################
#### Inputs
############################################################################################
results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/nakamura_hits'
hits_chchd2 <- read_excel("/Users/zacc/USyd/spatial_transcriptomics/analysis/nakamura_hits/hits_lists.xlsx", sheet = "CHCHD2_hits")
geomx_obj = "/Users/zacc/USyd/spatial_transcriptomics/analysis/nakamura_hits/geomx_110124_seurat.v5.Rdata" # geomx seurat object
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

############################################################################################
###### Part 1: Format Mouse DEGs
############################################################################################
setwd(results_folder)

# Format hits
#hits_chchd2$Gene <- toupper(hits_chchd2$Gene)
#hits_chchd2$Gene[hits_chchd2$Gene == "TRP53INP1"] <- "TP53INP1"
#hits_chchd2$Gene[hits_chchd2$Gene == "ZFP780B"] <- "ZNF780B"
#hits_chchd2 <- hits_chchd2[hits_chchd2$Gene != "JPX",] # note; could not find Jpx as it is a lncRNA 

# function to convert mouse to human
convert_mouse_to_human <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
}

musGenes <- hits_chchd2$Gene
hits_chchd2$Gene <- convert_mouse_to_human(musGenes)
targets <- hits_chchd2[hits_chchd2$Region %in% c("SN/VTA","SN"),]
###

############################################################################################
#### Part 2: Evaluation directional change of mouse DEGs in human SNpc ventral tier; Nanostring GeoMx
############################################################################################
# load GeoMx Seurat object and extract expression values and metadata
load(geomx_obj)

## evaluate all data
rdata <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata"
load(rdata)

# select counts and meta
counts <- as.matrix(gxdat_s@assays$RNA$data)
cell_meta <- as.data.frame(gxdat_s@meta.data)

# select subset
index <- which(cell_meta$ROI == "SNV" & cell_meta$Diagnosis %in% c("CTR","ILBD") & cell_meta$segment == "TH" )
counts <- counts[,index]
cell_meta <- cell_meta[index,]
#cell_meta$Diagnosis[cell_meta$Diagnosis != "CTR"] <- "PD"

colnames(cell_meta) <- make.names(colnames(cell_meta))
DEG <- targets$Gene[targets$Gene %in% row.names(counts)]
counts <- counts[DEG,]

# subset DEGs overlapping with counts
targets <- targets[targets$Gene %in% DEG,]

# # format metadata
cell_meta$Diagnosis <- as.factor(cell_meta$Diagnosis)
cell_meta$Age <- as.numeric(cell_meta$Age)
cell_meta$PMD.hs <- as.numeric(cell_meta$PMD.hs)
cell_meta$DV200 <- as.numeric(cell_meta$DV200)
cell_meta$plate <- unlist(lapply(strsplit(cell_meta$SampleID,"-"),function(x) x[2]))

## Run Voom
options(na.action='na.omit')
design <- model.matrix(~ Diagnosis +  Age + PMD.hs + DV200 + Sex + plate + Brainbank_ID,
                         #Brainbank_ID_recode,
                       data=cell_meta,
                       drop = FALSE)

# make names
colnames(design) <- make.names(colnames(design))

# create expression df & targ with design row.names
cell_meta <- cell_meta[row.names(design),]

# fit design
v <- voom(counts,design)
vfit <- lmFit(v)

# Perform LIMMA contrasts
cont.matrix <- makeContrasts(A="DiagnosisILBD",levels=design) # PD-CTR
fit2 <- contrasts.fit(vfit, cont.matrix)
vfit2 <- eBayes(fit2)
options(digits=3)

# Select significant DEGs and assign to list
lm_res <- topTable(vfit2,coef = "A", number = Inf,sort.by="P")

## Plotting 
set.seed(123)

# get directional changes and p-values add to targets
targets$Direction_geomx <- plyr::mapvalues(targets$Gene, from = row.names(lm_res), to = lm_res$logFC)
targets$Direction_geomx[targets$Direction_geomx < 0] <- "down"
targets$Direction_geomx[targets$Direction_geomx != "down"] <- "up"
targets$P.Value <- as.numeric(plyr::mapvalues(targets$Gene, from = row.names(lm_res), to = lm_res$adj.P.Val))

# Chi-square test
original <- targets$Direction
validation <- targets$Direction_geomx
overlap <- sum(original == validation)
# contingency table
table <- matrix(c(overlap, length(original) - overlap,
                  length(validation) - overlap, overlap), 
                nrow = 2)

# Chi-square test
chisq.test(table)

# p-value annots
pvalue = targets$P.Value
is_sig = pvalue < 0.05
pch = rep("*", length(is_sig))
pch[!is_sig] = NA
pvalue_col_fun = colorRamp2(c(0, 1,2, 3), c("grey", "white","yellow","orange")) # color mapping for -log10(pvalue)
# now we generate two legends, one for the p-value
lgd_pvalue = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3), 
                    labels = c("1", "0.1", "0.01", "0.001"))
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")

# row annots
row_ha = rowAnnotation(col = list(mouse=c("down" = "red", "up" = "green"),
                                  human=c("down" = "red", "up" = "green")),
                       mouse = targets$Direction,
                       human = targets$Direction_geomx,
                       pvalue = anno_simple(-log10(pvalue), col = pvalue_col_fun, pch = pch),
                       show_legend = FALSE)

# column annots
Subject_vec = viridis(length(unique(cell_meta$Brainbank_ID_recode)))
names(Subject_vec) <- unique(cell_meta$Brainbank_ID_recode)

age_fun = colorRamp2(range(cell_meta$Age), c("white", "red"))
dv200_fun = colorRamp2(range(cell_meta$DV200), c("white", "darkblue"))
pmd_fun = colorRamp2(range(cell_meta$PMD.hs), c("white", "darkgreen"))


col_ha = HeatmapAnnotation(
  Dx = cell_meta$Diagnosis,
  Subject = cell_meta$Brainbank_ID_recode,
  Sex = cell_meta$Sex,
  Age = cell_meta$Age,
  DV200 = cell_meta$DV200,
  PMD = cell_meta$PMD.hs,
  col = list(Dx = c("ILBD" = "hotpink1","CTR" = "black"),
             Subject = Subject_vec,
             Sex = c("M" = "lightblue","F" = "purple"),
             Age = age_fun,
             DV200 = dv200_fun,
             PMD = pmd_fun )
)

# plot heatmap
ht <- Heatmap(log10(counts), name = "Log10(Counts)", 
              right_annotation = row_ha,
              top_annotation = col_ha, 
              show_column_names = FALSE,
              cluster_columns = TRUE,
              row_names_gp = gpar(fontsize = 6))

# assign plot data
#name_plot = "hits_chchd2.SN.TH_n8"
#name_plot = "hits_chchd2.SN.Full_n8"
#name_plot = "hits_chchd2.VTA.TH_n7"
#name_plot = "hits_chchd2.VTA.Full_n7"
#name_plot = "hits_chchd2.SN.TH_ILBD.ePD.lPD_n8"
name_plot = "hits_chchd2.SN.ILBD.Full_n8"

# plot
pdf(paste0("heatmap_",name_plot,".pdf"))
draw(ht, annotation_legend_list = list(lgd_pvalue, lgd_sig))
dev.off()


############################################################################################
#### Part 3: tabulation of cohort statistics
############################################################################################
# create a glm of cohort stats
dplot <- cell_meta

dplot$Sex[dplot$Sex %in% c("male","Male")] <- "M"
dplot$Sex[dplot$Sex %in% c("female","Female")] <- "F"
dplot$Diagnosis <- as.factor(dplot$Diagnosis)
dplot$area <- as.numeric(dplot$area)
dplot$DV200 <- as.numeric(dplot$DV200)
dplot$Age <- as.numeric(dplot$Age)

# cases
tmp <- dplot

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
  dplyr::summarize(mean_age = median(Age,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(mean_age = IQR(Age,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(sd_age = sd(Age,na.rm=TRUE))

tmp2 <- tmp[!duplicated(tmp$Brainbank_ID),]
t.test(tmp2$Age[tmp2$Diagnosis == "CTR"],
       tmp2$Age[tmp2$Diagnosis == "ePD"])

#Sex
tmp %>%
  group_by(Diagnosis,Sex) %>%
  dplyr::summarise(count = n_distinct(Brainbank_ID))

tmp2 <- tmp[!duplicated(tmp$Brainbank_ID),]
chisq.test(tmp2$Sex,tmp2$Diagnosis)

# Dv200
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(result = median(DV200,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(result = IQR(DV200,na.rm=TRUE))

tmp2 <- tmp[!duplicated(tmp$Brainbank_ID),]
t.test(as.numeric(tmp2$DV200[tmp2$Diagnosis == "CTR"]),
       as.numeric(tmp2$DV200[tmp2$Diagnosis == "ePD"]))

# Area
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(mean_area = mean(area,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(sd_dv200 = sd(DV200,na.rm=TRUE))

# postmortem interval
tmp$`PMD hs` <- as.numeric(tmp$`PMD hs`)
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(result = median(`PMD hs`,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(result = IQR(`PMD hs`,na.rm=TRUE))


tmp2 <- tmp[!duplicated(tmp$Brainbank_ID),]
t.test(as.numeric(tmp2$`PMD hs`[tmp2$Diagnosis == "CTR"]),
       as.numeric(tmp2$`PMD hs`[tmp2$Diagnosis == "ePD"]))


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


############################################################################################
#### Part 4: Plot genes of interest
############################################################################################
genes_interest <- c("S100A6","SGK1","CALM2","SYT1","SNAP25","AP2A1","SYNGR3","TH","SLC6A3","SLC18A2")
genes_interest <- c("HSPA4","EIF3G")

# load geomx
load(geomx_obj)

# select counts and meta
counts <- as.matrix(gxdat_s@assays$RNA$data)
cell_meta <- as.data.frame(gxdat_s@meta.data)
colnames(cell_meta) <- make.names(colnames(cell_meta))
counts_2 <- counts[genes_interest,]
data_table <- cbind(cell_meta,t(counts_2))

res <- list()
y_list <- genes_interest
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "Diagnosis"
  y_variable = y_list[i]
  x_lab = "Diagnosis"
  y_lab = "Counts"
  colour_palette = c("ePD" = "red", "CTR" = "dodgerblue")
  
  # plot
  res[[i]] <- ggboxplot(
    data_table, x = x_variable, y = y_variable,
    fill = x_variable , palette = colour_palette, add = "dotplot") +
    ylab(y_lab) + xlab(x_lab) + theme(legend.position = "none") + ggtitle(y_variable) 
}

arrange <- ggarrange(plotlist=res, nrow=2, ncol=5, widths = c(2,2))
#ggsave("Box_genes_interest_140224.pdf", arrange)
#ggsave("Box_genes_interest_140224_2.pdf", arrange)


