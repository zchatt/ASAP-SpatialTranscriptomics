# Details: Scripts used for human validation of DEGs identified by chemogenetic (DREADD) mouse model to chronically hyperactivate of DA neurons.
# As described in manuscript "Chronic hyperactivation of midbrain dopamine neurons causes preferential dopamine neuron degeneration"

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
library(EnhancedVolcano)

############################################################################################
#### Inputs
############################################################################################
results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/nakamura_hits'
mouse_degs1 = "/Users/zacc/USyd/spatial_transcriptomics/analysis/nakamura_hits/hits_lists.xlsx" # interesting genes (n=14 + 2) + TH and SLC18A2
mouse_degs2 = "/Users/zacc/USyd/spatial_transcriptomics/analysis/nakamura_hits/240112_AllSNHits.xlsx" # 59 DEGs
geomx_obj = "/Users/zacc/USyd/spatial_transcriptomics/analysis/nakamura_hits/geomx_110124_seurat.v5.Rdata" # geomx seurat object
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

############################################################################################
###### Part 1: Format Mouse DEGs
############################################################################################
setwd(results_folder)

# Read in hits from AAV-Gq treated mice
# These are DEGs ID in both GqVeh (vehicle control) vs GqCNO (Gq receptor activation in neurons) and CNO only (activator control) vs GqCNO

## format gene list 1; interesting genes
hits_aavgq <- as.data.frame(read_excel(mouse_degs1, sheet = "AAV-Gq_hits"))
hits_aavgq$Gene <- row.names(hits_aavgq)
hits_aavgq$Gene <- toupper(hits_aavgq$Gene)
hits_aavgq$Gene[hits_aavgq$Gene == "DAT"] <- "SLC6A3"
hits_aavgq$Gene[hits_aavgq$Gene == "HBA-A2"] <- "HBA2"
hits_aavgq$Gene[hits_aavgq$Gene == "HBB-BT"] <- "HBA2"
hits_aavgq <- hits_aavgq[!duplicated(hits_aavgq$Gene),]
targets <- hits_aavgq[hits_aavgq$Region %in% c("VTA/SN","SN"),]

targets1 <- targets
targets1 <- rbind(targets,c("SN","down","SNCA"))
###

### 10 targets of interest 
targets <- targets[targets$Gene %in% c("SYT1","SLC6A3","TH","SYNGR3","CALM2","S100A6","SLC18A2","SLC10A4","SGK1"), ]
#targets <- rbind(targets,c("SN","down","SNCA"))
targets_alt <- targets$Gene
targets_alt[targets_alt == "SLC6A3"] <- "DAT"
targets_alt[targets_alt == "SLC18A2"] <- "VMAT2"


## format gene list 2; DEG 59
hits_aavgq <- as.data.frame(read_excel(mouse_degs2, sheet = "Sheet1"))
#hits_aavgq$Gene <- row.names(hits_aavgq)
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
musGenes <- hits_aavgq$Gene
hits_aavgq$Gene <- convert_mouse_to_human(musGenes)
targets <- hits_aavgq[hits_aavgq$Region %in% c("VTA/SN","SN"),]
targets2 <- targets
###


## only targets within dopamine and CA
#dopamine <- c("SYNGR3","SLC10A4","SLC18A2","TH","SLC6A3")
dopamine <- c("SLC10A4","SLC18A2")
#ca_transmission <- c("TRPC6","SYT1","CALM2","S100A6","SGK1")
ca_transmission <- c("SYT1","CALM2","SYNGR3")
targets$Group[targets$Gene %in% dopamine] <- "DA"
targets$Group[targets$Gene %in% ca_transmission] <- "Ca2+/Trans"
table(dopamine %in% targets$Gene)
table(ca_transmission %in% targets$Gene)
targets <- targets[!is.na(targets$Group),]

# get all genes for providing normalised data
targets3 <- rbind(targets1,targets2)
targets3 <- targets3[!duplicated(targets3$Gene),]


############################################################################################
#### Part 2: Evaluation directional change of mouse DEGs in human SNpc ventral tier; Nanostring GeoMx
############################################################################################
# load GeoMx Seurat object and extract expression values and metadata
load(geomx_obj)

# select counts and meta
counts <- as.matrix(gxdat_s@assays$RNA$data)
cell_meta <- as.data.frame(gxdat_s@meta.data)
colnames(cell_meta) <- make.names(colnames(cell_meta))
DEG <- targets$Gene[targets$Gene %in% row.names(counts)]
counts <- counts[DEG,]

# # or get genes of interest norm data
# counts <- as.matrix(gxdat_s@assays$RNA$data)
# DEG <- targets3$Gene[targets3$Gene %in% row.names(counts)]
# counts <- counts[DEG,]
# colnames(counts) == row.names(cell_meta)
# tmp <- cbind(cell_meta,t(counts))
# write.xlsx(tmp, file="human_SN.THpos_CTR.ePD_GeoMx_genesinterest_260224.xlsx")


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
design <- model.matrix(~ Diagnosis +  Age + PMD.hs + DV200 + Sex + plate + Brainbank_ID_recode,
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
cont.matrix <- makeContrasts(A="DiagnosisePD",levels=design) # PD-CTR
fit2 <- contrasts.fit(vfit, cont.matrix)
vfit2 <- eBayes(fit2)
options(digits=3)

# Select significant DEGs and assign to list
lm_res <- topTable(vfit2,coef = "A", number = Inf,sort.by="P")

# # volcano plot highlighting specific genes
# keyvals <- ifelse(
#   row.names(lm_res) %in% dopamine, 'purple',
#   ifelse(row.names(lm_res) %in% ca_transmission, 'royalblue',
#          'grey'))
# keyvals[is.na(keyvals)] <- 'black'
# names(keyvals)[keyvals == 'purple'] <- 'dopamine'
# names(keyvals)[keyvals == 'royalblue'] <- 'Ca2+ & transmission'
# 
# g1 <- EnhancedVolcano(lm_res,
#                       lab = rownames(lm_res),
#                       x = 'logFC',
#                       y = 'P.Value',
#                       selectLab = c(dopamine,ca_transmission),
#                       xlab = bquote(~Log[2]~ 'fold change'),
#                       pCutoff = NA,
#                       FCcutoff = NA,
#                       pointSize = 5,
#                       labSize = 4.5,
#                       # shape = c(6, 4, 2, 11),
#                       colCustom = keyvals,
#                       colAlpha = 1,
#                       legendPosition = 'left',
#                       legendLabSize = 15,
#                       legendIconSize = 5.0,
#                       drawConnectors = TRUE,
#                       widthConnectors = 1.0,
#                       colConnectors = 'black',
#                       arrowheads = FALSE,
#                       gridlines.major = TRUE,
#                       gridlines.minor = FALSE,
#                       border = 'partial',
#                       borderWidth = 1.5,
#                       borderColour = 'black') + ylim(0,1.5)
# 
# # save
# pdf(paste0("volcano_genes_interest_n10.pdf"))
# g1
# dev.off()

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
# row_ha = rowAnnotation(col = list(mouse=c("down" = "blue", "up" = "red"),
#                                   human=c("down" = "blue", "up" = "red")),
#                        mouse = targets$Direction,
#                        human = targets$Direction_geomx,
#                        pvalue = anno_simple(-log10(pvalue), col = pvalue_col_fun, pch = pch),
#                        show_legend = FALSE)

row_ha = rowAnnotation(col = list(mouse=c("down" = "blue", "up" = "red")),
                       mouse = targets$Direction,
                       show_legend = FALSE, 
                       simple_anno_size = unit(1, "cm"))
# column annots
Subject_vec = viridis(length(unique(cell_meta$Brainbank_ID_recode)))
names(Subject_vec) <- unique(cell_meta$Brainbank_ID_recode)

age_fun = colorRamp2(range(cell_meta$Age), c("white", "red"))
dv200_fun = colorRamp2(range(cell_meta$DV200), c("white", "darkblue"))
pmd_fun = colorRamp2(range(cell_meta$PMD.hs), c("white", "darkgreen"))


# col_ha = HeatmapAnnotation(
#   Dx = cell_meta$Diagnosis,
#   Subject = cell_meta$Brainbank_ID_recode,
#   Sex = cell_meta$Sex,
#   Age = cell_meta$Age,
#   DV200 = cell_meta$DV200,
#   PMD = cell_meta$PMD.hs,
#   col = list(Dx = c("ePD" = "hotpink1","CTR" = "black"),
#              Subject = Subject_vec,
#              Sex = c("M" = "lightblue","F" = "purple"),
#              Age = age_fun,
#              DV200 = dv200_fun,
#              PMD = pmd_fun )
# )

col_ha = HeatmapAnnotation(
  Dx = cell_meta$Diagnosis,
  col = list(Dx = c("ePD" = "hotpink1","CTR" = "black"))
)

# plot heatmap
ht <- Heatmap(log10(counts), name = "Log10(Counts)", 
              right_annotation = row_ha,
              top_annotation = col_ha, 
              show_column_names = FALSE,
              cluster_columns = TRUE,
              row_names_gp = gpar(fontsize = 6),
              row_split = targets$Group,
              cluster_row_slices = FALSE,
              row_title_gp = gpar(col = c("Dopamine" = "purple","Ca2+ & Transmission" = "orange"), cex = 0.7),
              width = ncol(counts)*unit(1, "mm"), 
              height = nrow(counts)*unit(5, "mm"))

# assign plot data
name_plot = "hits_aavgq_genes.interest"
#name_plot = "hits_aavgq.SN.TH_n59"

# plot
pdf(paste0("heatmap_",name_plot,".pdf"))
draw(ht, annotation_legend_list = list(lgd_pvalue, lgd_sig))
dev.off()

png(paste0("heatmap_",name_plot,".png"))
draw(ht, annotation_legend_list = list(lgd_pvalue, lgd_sig))
dev.off()

tiff(paste0("heatmap_",name_plot,".tiff"))
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


