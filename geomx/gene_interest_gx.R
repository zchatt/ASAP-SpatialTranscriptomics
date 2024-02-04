### Analysis of genes or gene-sets of interest from Spatial Transcriptomics (ST) and snRNAseq datasets

library(dplyr)
library(viridis)
library(data.table)
library(ggpubr)
library(rstatix)
library(spacexr)
library(Seurat)

###############
### Input
###############
# analysis dir
analysis_dir = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
setwd(analysis_dir)

# ST data
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata")

# ST cell deconvoluted data
myRCTD <- readRDS("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results/myRCTDde_gx.rds")

# DEG data
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results/cside_res_snvsnd.ctrpd.rds")

# single-cell data 
cell_ranger_data = "/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger" # single-cell cell ranger output
filenames <- list.files("/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger", pattern="*UMAP.tsv", full.names=TRUE) # read in metadata / umap of cell IDs
cell_id <- do.call(rbind, lapply(filenames,fread))
cell_id <- cell_id[!cell_id$NAME == "TYPE",]

# Genesets of interest
channel_members <- read_excel("/Users/zacc/USyd/spatial_transcriptomics/analysis/asap_edwards/channel_members_genelist.xlsx")
channel_members$Gene <- toupper(channel_members$Gene)

# source plotting functions
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/plotting_functions.R")

# colour palettes
Diagnosis_col = c("CTR"= "grey","ILBD" = "#00AFBB", "ePD" = "#E7B800","lPD" = "red")
age_group_col = magma(4)
ROI_col = c("SNV" = viridis(6)[1],
            "SNM" = viridis(6)[2],
            "SND" = viridis(6)[3],
            "SNL" = viridis(6)[4],
            "VTA" = viridis(6)[5],
            "LC" = viridis(6)[6])

############################################################################################
###### Part 1. Evaluating the expression of channel members in scRNAseq data
############################################################################################
# load single-cell data
sc_obj <- Read10X(cell_ranger_data) # single-cell cell ranger output

# create seurat and 
seurat_obj <- CreateSeuratObject(sc_obj)
seurat_obj <- SCTransform(seurat_obj)

# subset for cell-types
sc_obj <- sc_obj[,cell_id$NAME]
cell_types <- setNames(cell_id$Cell_Type, cell_id$NAME)

# get channel members
row.names(sc_obj)[row.names(sc_obj) %in% channel_members$Gene]
tmp <- sc_obj[row.names(sc_obj) %in% channel_members$Gene,]

table(cell_id$NAME == colnames(sc_obj))
cell_type <- as.data.frame(cell_id$Cell_Type)
data_table <- cbind(cell_type,t(tmp))
  
# create plot data
dplot <- as.data.frame(data_table) %>%
  group_by(`cell_id$Cell_Type`) %>% 
  summarise_all("median")
dplot1 <- as.matrix(t(dplot[,2:ncol(dplot)]))

# column annotations
Cell_Broad = gsub("_.*","",dplot$`cell_id$Cell_Type`)
anno_df = data.frame(
  Cell_Broad = Cell_Broad,
  Cell_type = dplot$`cell_id$Cell_Type`
)

ha = HeatmapAnnotation(df = anno_df)

# draw heatmap
pdf("heatmap_snrnaseq_channelmembers_median_annot.pdf", width=15, height=15)
hm <- Heatmap(dplot1, cluster_columns = FALSE,
              top_annotation = ha,
              column_split=as.factor(Cell_Broad)
              )
draw(hm,
     column_title = "Channel Members in Kamath et al",
     column_title_gp=grid::gpar(fontsize=16))
dev.off()


############################################################################################
###### Part 1b. Evaluating the expression of SNCA in scRNAseq & ST datasets
############################################################################################
## ST - GeoMx data ##
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat <- as.data.frame(gxdat_s@meta.data)

# i) Violin plots between ROIs for CTR 
gene_name <- "SNCA"
gene <- unlist(exp_dat[gene_name,])
dplot <- cbind(meta_dat,gene)
dplot <- dplot[dplot$segment == "TH" &  dplot$Diagnosis == "CTR",]

# calculate fold-changes
mean(dplot$gene[dplot$ROI == "SNV"]) / mean(dplot$gene[dplot$ROI == "SND"])
mean(dplot$gene[dplot$ROI == "SNV"]) / mean(dplot$gene[dplot$ROI == "LC"])
mean(dplot$gene[dplot$ROI == "SNV"]) / mean(dplot$gene[dplot$ROI == "VTA"])

# fact_lvls <- aggregate(gene ~ ROI, dplot, mean)$ROI[order(-aggregate(gene ~ ROI, dplot, mean)$gene)]
# dplot$ROI <- factor(dplot$ROI, levels = fact_lvls)
data_table <- dplot
colnames(data_table) <- make.names(colnames(data_table))

dplot2 <- dplot[,c("ROI","gene")]
v1 <- violin_plot_function(dplot2,"ROI","gene","ROI",paste0(gene_name, " expression"), ROI_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="CTR ROI's (TH)")


# linear mixed-effects model of SNCA
y_list <- c("gene")
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "ROI"
  y_variable = y_list[i]
  x_lab = "ROI"
  y_lab = "Normalized counts"
  colour_palette = ROI_col 
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID) + (1 | Age) + (1 | Sex) + (1 | PMD.hs)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(ROI = "Tukey"))
  summary(posthoc)
  
  # format for plotting
  tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  
  # make violin plot 
  dplot2 <- data_table[,c(x_variable,y_variable)]
  bxp <- ggviolin(
    dplot2, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab(y_lab) + xlab(x_lab) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                     geom="pointrange", color="black")
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  
  # get sig values
  if(length(stat.test$p.adj < 0.05) > 1){
    # add p-values to plot and save the result
    stat.test <- stat.test[stat.test$p.adj < 0.05, ]
    bxp <- bxp + stat_pvalue_manual(stat.test,
                                    y.position = max(data_table$y) * 1.6, step.increase = 0.1,
                                    label = "p.adj.signif") 
  }
  
  # save plot
  arrange <- ggarrange(plotlist=list(bxp), nrow=2, ncol=2)
  ggsave(paste("violin_SNCA.CTR.ROI_nothresh",y_variable,".png"), arrange)
  
}


#i.b) violin plot with iNM and SNCA
# format iNM
quant_data = "/Users/zacc/USyd/NM_analysis/df_iNMeNM_230124.Rdata"
meta_data = "/Users/zacc/USyd/NM_analysis/NM_data_170124/ASAP1-IHC cohort.xlsx"
load(quant_data)
meta <- read_xlsx(meta_data,1)
meta$Diagnosis[meta$Diagnosis == "Ct"] <- "CTR"
df_agg_merged <- merge(df_agg,meta,by="Brainbank_ID",all.x=TRUE)
# evaluate counts/ um2
data_table <- df_agg_merged[df_agg_merged$intra.extra == "iNM" ,]
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID_recode, ROI,Age,Sex,PMD,`ROI Area [µm²]`, sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$n_per.µm2 <- dplot$n / as.numeric(dplot$`ROI Area [µm²]`)
d1 <- dplot[complete.cases(dplot$ROI),c("Diagnosis","ROI","n_per.µm2")]
colnames(d1) <- c("Diagnosis","ROI","value")
d1$Measurement <- c("iNM (n/µm²)")

mean_ctr <- mean(d1$value[d1$Diagnosis == "CTR"])
d1$value <- d1$value/mean_ctr * 100

# format SNCA
gene_name <- "SNCA"
gene <- exp_dat[gene_name,]
dplot <- cbind(meta_dat,gene)
d2 <- dplot[dplot$segment == "TH",c("Diagnosis","ROI","gene")]
colnames(d2) <- c("Diagnosis","ROI","value")
d2$Measurement <- c("SNCA mRNA")

mean_ctr <- mean(d2$value[d2$Diagnosis == "CTR"])
d2$value <- d2$value/mean_ctr * 100

# combine and plot
dplot <- rbind(d1,d2)
dplot <- dplot[dplot$ROI == "SNV",]
dplot$Diagnosis <- factor(dplot$Diagnosis,levels=c("CTR","ILBD", "ePD","lPD"))

bxp <- ggviolin(dplot, x = "Diagnosis", y = "value", fill = "Measurement", palette = c("pink","grey")) + 
  geom_boxplot(aes(color = Measurement), width = 0.15, position = position_dodge(0.8)) +
  scale_color_manual(values = c("black","black")) + ylab("Change relative to CTR (%) ") + geom_hline(yintercept = 100, lty = 2) 


# save plot
arrange <- ggarrange(plotlist=list(bxp), nrow=2, ncol=2)
ggsave(paste("violin_joint.SNCA.iNM.SNV.png"), arrange)







# ii) high/ low SNCA expression cell-types
## ST - GeoMx data ##
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat <- as.data.frame(gxdat_s@meta.data)

# i) Violin plots between ROIs for CTR 
gene_name <- "SNCA"
gene <- exp_dat[gene_name,]
dplot <- cbind(meta_dat,gene)
#dplot2 <- dplot[dplot$ROI == "SNV" ,]
dplot2 <- dplot[dplot$Brainregion == "A9" & dplot$segment == "TH",]
#dplot2 <- dplot[dplot$Brainregion == "A9",]
dplot2$high_SNCA <- dplot2$gene > 250

# get cell-type proportions and subset with 
norm_weights <- normalize_weights(myRCTD@results$weights)
dplot3 <- merge(dplot2,norm_weights, by = "row.names")
dplot3$Diagnosis <- as.factor(dplot3$Diagnosis)
dplot3$segment <- as.factor(dplot3$segment)

# stats













# iii) plots of cells of interest
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat <- as.data.frame(gxdat_s@meta.data)
dplot2 <- meta_dat[meta_dat$segment == "TH" & meta_dat$ROI == "SNV" ,]

# get cell-type proportions and format
norm_weights <- normalize_weights(myRCTD@results$weights)
dplot3 <- merge(dplot2,norm_weights, by = "row.names")
dplot3$Diagnosis <- as.factor(dplot3$Diagnosis)
dplot3$segment <- as.factor(dplot3$segment)
dplot3$Diagnosis <- factor(dplot3$Diagnosis,levels=c("CTR","ILBD", "ePD","lPD"))
dplot3$Diagnosis_CTR <- as.factor(dplot3$Diagnosis == "CTR")
data_table <- dplot3
colnames(data_table) <- make.names(colnames(data_table))

data_table$MG_Inactive <- rowSums(data_table[,c("MG_TSPO_VIM","MG_CECR2_FGL1")])
data_table$MG_Active <- rowSums(data_table[,c("MG_CCL3","MG_GPNMB_LPL","MG_FOSL2","MG_GPNMB_SULT1C2","MG_GPNMB_SULT1C2","MG_GPNMB_SUSD1")])
data_table$MG_AI_ratio <- data_table$MG_Active / data_table$MG_Inactive # calculate a microglia activation


# NOTES:
## potentially inactive microglia subtypes
# MG_TSPO_VIM; TSPO is a marker of activated microglia in rodent BUT NOT human https://pubmed.ncbi.nlm.nih.gov/37640701/
# MG_CECR2_FGL1; CECR2 is lower expression in age and AD microglia - https://www.cell.com/cell-reports/pdfExtended/S2211-1247(20)30824-X & lower expression in LPS treated microglia https://www.nature.com/articles/s41398-022-02265-6

## potentially activated microglia subtypes
# MG_CCL3; CCL3 is a marker of activated microglia in humans - https://pubmed.ncbi.nlm.nih.gov/37358357/
# MG_GPNMB_LPL, MG_FOSL2, MG_GPNMB_SULT1C2, MG_GPNMB_SULT1C2, MG_GPNMB_SUSD1; GPNMB, LPL, FOSL2 is expressed along microglial activation trajectory in midbrain PD - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9050543/#sup1
# MG_SPON1; SPON1 is a marker of an inflammatory state of microglia id in CTR and AD of several brain regions - https://www.cell.com/cell/pdf/S0092-8674(23)00971-6.pdf

## cell cycling
# MKI67; mouse data shows increased following nerve injury - https://www.nature.com/articles/s41421-022-00377-3, but is a cell-cycling gene - https://www.sciencedirect.com/science/article/pii/S221112472030019X#bib33

## general expression in microglia
# MG_OPRM1; OPRM1 appears to be ubiquitously expressed in microglia - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6328486/
 
## unknwon
# "MG_MGAM"               

# load published gene sets, include microglia activation trajectories
load("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/deg_markers/tables/deg_list.published.R")
tmp <- deg_list[["Smajić.2022.Microglia_TG"]]

# linear mixed-effects model of RGB between Ctr and ILBD within each ROI; area
#y_list <- c("MG_CCL3","MG_MGAM","MG_TSPO_VIM")
#y_list <- colnames(data_table)[grep("MG_",colnames(data_table))]
y_list <- c("MG_MGAM","MG_AI_ratio")
res <- list()
roi_uniq <- unique(data_table$ROI)
for(z in 1:length(roi_uniq)){
  print(z)
  roi_cont <- roi_uniq[z]
  data_table2 <- data_table[data_table$ROI == roi_cont,]
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "Diagnosis"
  y_variable = y_list[i]
  x_lab = "Diagnosis"
  y_lab = "Relative to CTR (%)"
  colour_palette = Diagnosis_col
  
  # get % relative to control
  data_table2$value <- data_table2[,y_variable]
  mean_ctr <- mean(data_table2$value[data_table2$Diagnosis == "CTR"])
  data_table2$value <- data_table2$value/mean_ctr * 100
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste("value","~",x_variable," + (1 | Brainbank_ID) + (1 | Age) + (1 | Sex) + (1 | PMD.hs)")), data = data_table2)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
  summary(posthoc)
  
  # make violin plot 
  tmp <- data_table2[,c(eval(x_variable),"value")]
  bxp <- ggviolin(
    tmp, x = x_variable, y = "value", 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    ylab(y_lab) + xlab(x_lab) + stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="black") +
    geom_hline(yintercept = 100, lty = 2) 
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  # get sig values
  if(any(stat.test$p.adj < 0.05)){
    print("test")
    # add p-values to plot and save the result
    stat.test <- stat.test[stat.test$p.adj < 0.05, ]
    bxp <- bxp + stat_pvalue_manual(stat.test,
                                    y.position = max(data_table2$y) * 0.8, step.increase = 0.1,
                                    label = "p.adj.signif") 
  }
  
  
  # save plot
  arrange <- ggarrange(plotlist=list(bxp), nrow=3, ncol=2, widths = c(2,2))
  ggsave(paste("Vln_CellProp.A9.TH.",y_variable,".",x_lab,".",roi_cont,".png"), arrange)
}
}



















data_table <- gather(dplot3, condition, measurement, Astro_CYP4F12:SOX6_PART1, factor_key=TRUE)
data_table <- data_table[,c("high_SNCA","condition","measurement")]

stat.test <- data_table %>%
  group_by(condition) %>%
  t_test(measurement ~ high_SNCA) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

# split by p-value
stat.test1 <- stat.test[stat.test$p.adj < 0.05,]
stat.test2 <- stat.test[stat.test$p.adj > 0.05,]

data_table1 <- data_table[data_table$condition %in% stat.test1$condition,]
data_table2 <- data_table[data_table$condition %in% stat.test2$condition,]

# boxplot
grp_col <- c("grey","black")
b1 <- boxplot_2grp_function(data_table1,"condition","measurement","Cell-Type","Cell-Type Proportion","high_SNCA",grp_col) +
  theme(plot.title = element_text(hjust = 1)) +  rotate_x_text(45) +
  labs(title="CTR ROI's (TH): high v low SNCA")

arrange <- ggarrange(plotlist=list(b1), nrow=2, ncol=1)
ggsave("Box_CTR.TH.ROI.SNCA.highlow_BHp0.05_celltypes.png", arrange)

b1 <- boxplot_2grp_function(data_table2,"condition","measurement","Cell-Type","Cell-Type Proportion","high_SNCA",grp_col) +
  theme(plot.title = element_text(hjust = 1)) +  rotate_x_text(45) +
  labs(title="CTR ROI's (TH): high v low SNCA")

arrange <- ggarrange(plotlist=list(b1), nrow=2, ncol=1)
ggsave("Box_CTR.TH.ROI.SNCA.highlow_BHpnonsig_celltypes.png", arrange)






## Violin by diagnosis for cells of interest
gene <- exp_dat[gene_name,]
dplot <- cbind(meta_dat,gene)
dplot2 <- merge(dplot,norm_weights, by = "row.names")
#dplot2 <- dplot2[dplot2$gene > 250,]
dplot2 <- dplot2[dplot2$ROI == "SNV" &  dplot2$segment == "TH",]

dplot2$Diagnosis <- factor(dplot2$Diagnosis, levels = c("CTR","ILBD","ePD","lPD"))
data_table <- dplot2
colnames(data_table) <- make.names(colnames(data_table))

# linear mixed-effects model & violin
y_list <- c("MG_MGAM","MG_TSPO_VIM")
for(i in 1:length(y_list)){
  print(i)
  # define variables
  x_variable = "Diagnosis"
  y_variable = y_list[i]
  x_lab = "Diagnosis"
  y_lab = paste0("Microglia (",y_variable,") proportion")
  colour_palette = Diagnosis_col
  
  # Fit a linear mixed-effects model
  model <- lmer(formula(paste(eval(y_variable),"~",x_variable," + (1 | Brainbank_ID) + (1 | Age) + (1 | Sex) + (1 | PMD.hs)")), data = data_table)
  
  # perform post-hoc tests to compare different regions using Tukey method
  posthoc <- glht(model, linfct = mcp(Diagnosis = "Tukey"))
  summary(posthoc)
  
  # # format for plotting
  # tmp <- aggregate(formula(paste(y_variable," ~ ",x_variable)), data_table, mean)
  # fact_lvls <- tmp[order(-tmp[,y_variable]),][,x_variable]
  # data_table[,x_variable] <- factor(data_table[,x_variable], levels = fact_lvls)
  # 
  # make violin plot 
  dplot2 <- data_table[,c(x_variable,y_variable)]
  bxp <- ggviolin(
    dplot2, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylab(y_lab) + xlab(x_lab) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + stat_summary(fun.data=mean_sdl, mult=1, 
                                                                                     geom="pointrange", color="black")
  
  # set x and y in data_table for plotting
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  
  # create stat.test df
  summary_posthoc <- summary(posthoc)
  group1 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[1])))
  group2 = gsub(" ","",unlist(lapply(strsplit(names(summary_posthoc$test$coefficients),"-"), function(x) x[2])))
  p.adj = as.numeric(summary_posthoc$test$pvalues)
  stat.test <- as.data.frame(cbind(group1,group2,p.adj))
  stat.test$p.adj <- as.numeric(stat.test$p.adj)
  stat.test$p.adj.signif <- symnum(stat.test$p.adj, corr = FALSE, na = FALSE, 
                                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                   symbols = c("***", "**", "*", ".", " "))
  print(stat.test)
  
  # get sig values
  if(length(stat.test$p.adj < 0.05) > 1){
    # add p-values to plot and save the result
    stat.test <- stat.test[stat.test$p.adj < 0.05, ]
    bxp <- bxp + stat_pvalue_manual(stat.test,
                                    y.position = max(data_table$y) * 1.6, step.increase = 0.1,
                                    label = "p.adj.signif") 
  }
  
  # save plot
  arrange <- ggarrange(plotlist=list(bxp), nrow=2, ncol=2)
  ggsave(paste("violin_SNCA.CTR.ROI",y_variable,".png"), arrange)
  
}






bxp1 <- ggboxplot(
  dplot3,"Diagnosis","Ex_LAMP5_NTNG2", palette = Diagnosis_col) + 
  theme(legend.position = "none") + ylab("cell proportion") +
  rotate_x_text(65) + labs(title="Ex_LAMP5_NTNG2: SNV TH ROI's") +
  geom_pwc(
    aes(group = Diagnosis), tip.length = 0,
    method = "t_test", label = "p.adj.signif", p.adjust.method = "BH"
  )

bxp2 <- ggboxplot(
  dplot2,"Diagnosis","Ex_LAMP5_BAIAP3", palette = Diagnosis_col) + 
  theme(legend.position = "none") + ylab("cell proportion") +
  rotate_x_text(65) + labs(title="Ex_LAMP5_BAIAP3: SNV TH ROI's") +
  geom_pwc(
    aes(group = Diagnosis), tip.length = 0,
    method = "t_test", label = "p.adj.signif", p.adjust.method = "BH"
  )


bxp3 <- ggboxplot(
  dplot2,"Diagnosis","MG_MGAM", palette = Diagnosis_col) + 
  theme(legend.position = "none") + ylab("cell proportion") +
  rotate_x_text(65) + labs(title="MG_MGAM: SNV TH ROI's") +
  geom_pwc(
    aes(group = Diagnosis), tip.length = 0,
    method = "t_test", label = "p.adj.signif", p.adjust.method = "BH"
  )

bxp4 <- ggboxplot(
  dplot2,"Diagnosis","MG_TSPO_VIM", palette = Diagnosis_col) + 
  theme(legend.position = "none") + ylab("cell proportion") +
  rotate_x_text(65) + labs(title="MG_TSPO_VIM: SNV TH ROI's") +
  geom_pwc(
    aes(group = Diagnosis), tip.length = 0,
    method = "t_test", label = "p.adj.signif", p.adjust.method = "BH"
  )

arrange <- ggarrange(plotlist=list(bxp1,bxp4,bxp2,bxp3), nrow=2, ncol=4)
ggsave("box_celltypes_SNV.TH.png", arrange)

# iii) Violin plots between CTR & PD for each ROI
gene_name <- "SNCA"
gene <- exp_dat[gene_name,]
dplot <- cbind(meta_dat,gene)
dplot <- dplot[dplot$segment == "TH" & dplot$ROI != "SNL",]
dplot$Diagnosis <- factor(dplot$Diagnosis, levels = c("CTR","ILBD","ePD","lPD"))

# plot violin
roi <- unique(dplot$ROI)
res <- list()
for (i in 1:length(roi)){
  print(roi[i])
  dplot2 <- dplot[dplot$segment == "TH" & dplot$ROI == roi[i],c("Diagnosis","gene")]
  dplot2$Diagnosis <- droplevels(dplot2$Diagnosis)
  
  res[[i]] <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
    theme(plot.title = element_text(hjust = 1)) +
    labs(title=roi[i])
}

arrange <- ggarrange(plotlist=res[1:4], nrow=2, ncol=2)
ggsave("Violin_TH.ROI.SNCA.Dx.png", arrange)


## scRNAseq - Kamath data ##
# get gene interest
tmp <- as.data.frame(sc_obj[row.names(sc_obj) %in% c("SNCA","TH"),])

table(cell_id$NAME == colnames(sc_obj))
cell_type <- as.data.frame(cell_id$Cell_Type)
data_table <- cbind(cell_type,t(tmp))

# plot SCNA expression
dplot <- data_table
fact_lvls <- aggregate(SNCA ~ `cell_id$Cell_Type`, dplot, mean)$`cell_id$Cell_Type`[order(-aggregate(SNCA ~ `cell_id$Cell_Type`, dplot, mean)$SNCA)]
dplot$`cell_id$Cell_Type` <- factor(dplot$`cell_id$Cell_Type`, levels = fact_lvls)

# make boxplot
colnames(dplot) <- c("cell","SNCA","TH")
bxp <- ggboxplot(
  dplot,"cell","SNCA", palette = colour_palette, order = fact_lvls) + 
  theme(legend.position = "none") +
  rotate_x_text(65) + labs(title="Kamath et al snRNAseq")

arrange <- ggarrange(plotlist=list(bxp), nrow=2, ncol=1)
ggsave("boxplot_kamath_SCNA.png", arrange)


dplot2 <- dplot[dplot$cell %in% c("SOX6_DDT","Ex_LAMP5_BAIAP3","Ex_LAMP5_NTNG2","MG_MGAM","MG_CCL3","MG_TSPO_VIM"),]

dplot2$cell <- droplevels(dplot2$cell)
stat.test <- dplot2 %>% t_test(SNCA ~ cell) %>% add_significance("p.adj")
stat.test <- stat.test[stat.test$p.adj < 0.05, ]
print(stat.test)

# add p-values to plot and save the result
bxp <- bxp + stat_pvalue_manual(stat.test,
                                y.position = max(data_table$y) * 1.6, step.increase = 0.1,
                                label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)



## boxplots of cells of interest
bxp1 <- ggviolin(
  dplot2,"cell","SNCA", palette = colour_palette, order = fact_lvls) + 
  theme(legend.position = "none") +
  rotate_x_text(65) + labs(title="Kamath et al snRNAseq") +
  geom_pwc(
    aes(group = cell), tip.length = 0,
    method = "t_test", label = "p.adj.signif", p.adjust.method = "BH"
  )

dplot2 <- dplot[dplot$cell %in% c("SOX6_DDT","MG_MGAM","MG_CCL3","MG_TSPO_VIM"),]
bxp2 <- ggviolin(
  dplot2,"cell","SNCA", palette = colour_palette, order = fact_lvls) + 
  theme(legend.position = "none") +
  rotate_x_text(65) + labs(title="Kamath et al snRNAseq") +
  geom_pwc(
    aes(group = cell), tip.length = 0,
    method = "t_test", label = "p.adj.signif", p.adjust.method = "BH"
  )

arrange <- ggarrange(plotlist=list(bxp1,bxp2), nrow=1, ncol=4)
ggsave("boxplot_kamath_SCNA_cellselect.png", arrange)


############################################################################################
###### Part 1b. Evaluating the expression of NM biosynthesis members
############################################################################################
NM_synthesis_genes <- c("TH","TYR","DDC","COMT","MAOA","MAOB","AKR1B1","SLC18A2","RGN",row.names(sc_obj)[grep("ALDH",row.names(sc_obj))])

# get channel members
tmp <- sc_obj[row.names(sc_obj) %in% NM_synthesis_genes,]

table(cell_id$NAME == colnames(sc_obj))
cell_type <- as.data.frame(cell_id$Cell_Type)
data_table <- cbind(cell_type,t(tmp))

# create plot data
dplot <- as.data.frame(data_table) %>%
  group_by(`cell_id$Cell_Type`) %>% 
  summarise_all("median")
dplot1 <- as.matrix(t(dplot[,2:ncol(dplot)]))

# column annotations
Cell_Broad = gsub("_.*","",dplot$`cell_id$Cell_Type`)
anno_df = data.frame(
  Cell_Broad = Cell_Broad,
  Cell_type = dplot$`cell_id$Cell_Type`
)

ha = HeatmapAnnotation(df = anno_df)

# draw heatmap
pdf("heatmap_snrnaseq_NMsynthesisgenes.pdf", width=15, height=15)
hm <- Heatmap(dplot1, cluster_columns = FALSE,
              top_annotation = ha,
              column_split=as.factor(Cell_Broad)
)
draw(hm,
     column_title = "NM synthesis genes in Kamath et al",
     column_title_gp=grid::gpar(fontsize=16))
dev.off()


############################################################################################
###### Part 2. Printing gene expression of interest to excell
############################################################################################
## GeoMx data ##
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat <- as.data.frame(gxdat_s@meta.data)

# Get gene expression for Nakamura Lab
gene_name <- c("CHCHD10","CHCHD2","SNCA","PDHA1")
gene <- exp_dat[gene_name,]
dplot <- cbind(meta_dat,t(gene))
dplot <- as.data.frame(dplot)
dplot <- dplot[dplot$segment == "TH" & dplot$ROI %in% c("SNV","SND"),]
dplot <- dplot[,c("Brainbank_ID","Sex","Age","PMD hs","ROI","SNCA","PDHA1","CHCHD10","CHCHD2")]
openxlsx::write.xlsx(dplot, file = "targeted_genelist_geomx_011223.xlsx")

# Get gene expression for YF
gene_name <- c("KCNJ6", "SLC6A3", "SLC18A2", "TH", "ALDH1A1","SNCA")
gene <- exp_dat[gene_name,]
dplot <- cbind(meta_dat,t(gene))
dplot <- as.data.frame(dplot)
dplot <- dplot[dplot$segment == "TH" & dplot$ROI %in% c("SNV","SND","VTA"),]
dplot <- dplot[,c(c("Brainbank_ID","Sex","Age","PMD hs","ROI"),gene_name)]
openxlsx::write.xlsx(dplot, file = "targeted_genelist_geomx_151223.xlsx")



############################################################################################
###### Part 3. Violin plots of PGK1
############################################################################################
# geomx data
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat <- as.data.frame(gxdat_s@meta.data)
gene_name <- "PGK1"
gene <- exp_dat[gene_name,]
dplot <- cbind(meta_dat,gene)
dplot$Diagnosis <- factor(dplot$Diagnosis, levels = c("CTR","ILBD","ePD","lPD"))

# plot all roi's independently
dplot <- dplot[dplot$segment == "Full ROI" & dplot$ROI != "RN",]
roi <- unique(dplot$ROI)
res <- list()
for (i in 1:length(roi)){
  print(roi[i])
  dplot2 <- dplot[dplot$ROI == roi[i],c("Diagnosis","gene")]
  dplot2$Diagnosis <- droplevels(dplot2$Diagnosis)
  
  res[[i]] <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
    theme(plot.title = element_text(hjust = 1)) +
    labs(title=roi[i])
}

arrange <- ggarrange(plotlist=res[1:4], nrow=2, ncol=2)
ggsave("Violin_Full.ROI_Dx_PGK1.png", arrange)
arrange <- ggarrange(plotlist=res[5:6], nrow=2, ncol=2)
ggsave("Violin_Full.ROI_Dx_PGK1_2.png", arrange)





















############################################################################################
###### Part X: Post C-SIDE; GSEA and plotting
############################################################################################

# load CSIDE results; generated from ctDEG_gx, list of DEG's should be named "cside_res"
load("cside_res_snvsnd.ctrpd.rds")

# define gene sets
# NOTE: using "H: hallmark gene sets", using "C2: curated gene sets", using "C5: ontology gene sets", using "C8: cell type signature gene sets"
all_gene_sets = msigdbr(species = "human")
msigdbr_df <- all_gene_sets %>%
  dplyr::filter(gs_cat == "C2")
msigdbr_list_c2 = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
msigdbr_df <- all_gene_sets %>%
  dplyr::filter(gs_cat == "C5")
msigdbr_list_c5 = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
msigdbr_df<- all_gene_sets %>%
  dplyr::filter(gs_cat == "C8")
msigdbr_list_c8 = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# run per DEG list
res_tmp <- list()
for (i in 1:length(cside_res)){
  print(i)
  
  # reduce to sig DEG's
  sig_gene_list <- cside_res[[i]]
  sig_gene_list <- lapply(sig_gene_list, function(x) x[x$p_val < 0.05,])
  sig_gene_list <- lapply(sig_gene_list, function(x) x[order(x$p_val),])
  res_tmp[[i]] <- unlist(lapply(sig_gene_list, row.names))
  
  ## Plots
  # barplot for n sig genes x cell-type
  x = unlist(lapply(sig_gene_list, nrow))
  
  dplot <- as.data.frame(cbind(cell_type=names(x),sig_genes=x))
  dplot$sig_genes <- as.numeric(dplot$sig_genes)
  
  g1 <- ggplot(dplot, aes(x=reorder(cell_type,-x),y=x)) +
    geom_bar(stat="identity", width=0.7)+
    theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
    xlab("Cell Type") + ylab("Sig CtDEG's (n)")
  
  # save plots
  ggsave(paste0("barplot_ctdeg_",name1,".png"),
         ggarrange(g1,nrow=2,ncol=2), device = "png")
  
  # all DEGs for FGEA
  sig_gene_list <- cside_res[[i]]
  sig_gene_list <- lapply(sig_gene_list, function(x) x[order(x$p_val),])
  
  for (z in 1:length(sig_gene_list)){
    # select res from list, perform FGSEA and plot results
    print(names(sig_gene_list)[z])
    res <- as.data.frame(sig_gene_list[z])
    name1 <- make.names(names(sig_gene_list)[z])
    gene_stats <- res[,grep("log_fc",colnames(res))]
    names(gene_stats) <- row.names(res)
    
    # FGSEA - C2
    fgseaRes <- fgsea(pathways = msigdbr_list_c2,
                      stats    = gene_stats,
                      minSize  = 15,
                      maxSize  = 500)
    fgseaRes <- fgseaRes[order(pval), ]
    head(fgseaRes)
    
    # plot
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    pdf(paste0("FGSEA_C2_",name1,".pdf"))
    print(plotGseaTable(msigdbr_list_c2[topPathways], gene_stats, fgseaRes,
                        gseaParam=0.5, pathwayLabelStyle=list(size=6, color="black")))
    dev.off()
    
    # FGSEA - C5
    fgseaRes <- fgsea(pathways = msigdbr_list_c5,
                      stats    = gene_stats,
                      minSize  = 15,
                      maxSize  = 500)
    fgseaRes <- fgseaRes[order(pval), ]
    head(fgseaRes)
    
    # plot
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    pdf(paste0("FGSEA_C5_",name1,".pdf"))
    print(plotGseaTable(msigdbr_list_c5[topPathways], gene_stats, fgseaRes,
                        gseaParam=0.5, pathwayLabelStyle=list(size=6, color="black")))
    dev.off()
    
    # FGSEA - C8
    fgseaRes <- fgsea(pathways = msigdbr_list_c8,
                      stats    = gene_stats,
                      minSize  = 15,
                      maxSize  = 500)
    fgseaRes <- fgseaRes[order(pval), ]
    head(fgseaRes)
    
    # plot
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    pdf(paste0("FGSEA_C8_",name1,".pdf"))
    print(plotGseaTable(msigdbr_list_c8[topPathways], gene_stats, fgseaRes,
                        gseaParam=0.5, pathwayLabelStyle=list(size=6, color="black")))
    dev.off()
  }
}









############################################################################################
###### Part 3: Examing DaN subtypes with unique function
############################################################################################

Aldh1a1_pos <- c(toupper(c("Aldh1a1"))) 
Aldh1a1_Anxa1_pos <- c(toupper(c("Aldh1a1","Anxa1"))) #innervate a more dorsally restricted region of the striatum  than Aldha1a1 alone
Vglut2_pos <- c(toupper(c("Vglut2","Calb1+")))
Calb1_pos <- c(toupper(c("Calb1")))

Aldh1a1 <- c(toupper(c()))                 
Aldh1a1_Anxa1_neg <- c(toupper(c()))            
Vglut2_neg <- c(toupper(c("Otx2","Sox6","Aldh1a1","Crhbp")))
Calb1_neg <- c(toupper(c("Vglut2","Otx2","Sox6","Aldh1a1")))
                


# plot the kamath single-cell data to identify analagous neuron subtypes than identified in Azcorra et al.
DotPlot(pbmc3k.final, features = features, split.by = "groups") + RotatedAxis()


             
### mouse dopamine neuron subtypes
genes_interest <- toupper(c("Aldh1a1","Slc17a6","Calb1","Otx2","Sox6","Calb1","Crhbp","TH"))


#Vglut2+ = Vglut2+/Otx2−/Sox6−/Aldh1a1−/Calb1+/Crhbp−,
#Calb1+ = Calb1+/Vglut2−/Otx2−/Sox6−/Aldh1a1−.

# get channel members
tmp <- sc_obj[row.names(sc_obj) %in% genes_interest,]

table(cell_id$NAME == colnames(sc_obj))
cell_type <- as.data.frame(cell_id$Cell_Type)
data_table <- cbind(cell_type,t(tmp))

# create plot data
dplot <- as.data.frame(data_table) %>%
  group_by(`cell_id$Cell_Type`) %>% 
  summarise_all("median")
dplot1 <- as.matrix(t(dplot[,2:ncol(dplot)]))

# column annotations
Cell_Broad = gsub("_.*","",dplot$`cell_id$Cell_Type`)
anno_df = data.frame(
  Cell_Broad = Cell_Broad,
  Cell_type = dplot$`cell_id$Cell_Type`
)

ha = HeatmapAnnotation(df = anno_df)

# draw heatmap
#pdf("heatmap_snrnaseq_channelmembers_median_annot.pdf", width=15, height=15)
hm <- Heatmap(dplot1, cluster_columns = FALSE,
              top_annotation = ha,
              column_split=as.factor(Cell_Broad)
)
draw(hm,
     column_title = "Genes interest in Kamath et al",
     column_title_gp=grid::gpar(fontsize=16))
#dev.off()


