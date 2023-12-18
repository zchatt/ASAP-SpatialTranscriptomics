### Analysis of genes of interest from GeoMx and snRNAseq datasets

library(dplyr)
library(viridis)
library(data.table)
library(ggpubr)
library(rstatix)
library(spacexr)
library(Seurat)


### Inputs
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
setwd(analysis_dir)

# ST data
rdata = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata"

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
###### Part 1b. Evaluating the expression of SNCA in scRNAseq & ST
############################################################################################

## ST ##
# colour palettes
Diagnosis_col = c("grey","#00AFBB", "#E7B800","red")
ROI_col = viridis(6)

# Violin plot function
violin_plot_function <- function(data_table, x_variable, y_variable, x_lab, y_lab, colour_palette) {
  # make violin plot
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none")
  
  # perform t-test
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  stat.test <- data_table %>% t_test(y ~ x) %>% add_significance("p.adj")
  stat.test <- stat.test[stat.test$p.adj < 0.05, ]
  print(stat.test)
  
  # add p-values to plot and save the result
  bxp <- bxp + stat_pvalue_manual(stat.test,
                                  y.position = max(data_table$y) * 1.6, step.increase = 0.1,
                                  label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
  
  #return object
  return(bxp)
}

## GeoMx data ##
load(rdata)
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat <- as.data.frame(gxdat_s@meta.data)

# i) Violin plots between ROIs for CTR 
gene_name <- "SNCA"
gene <- exp_dat[gene_name,]
dplot <- cbind(meta_dat,gene)
dplot <- dplot[dplot$segment == "TH" &  dplot$Diagnosis == "CTR",]

# calculate fold-changes
mean(dplot$gene[dplot$ROI == "SNV"]) / mean(dplot$gene[dplot$ROI == "SND"])
mean(dplot$gene[dplot$ROI == "SNV"]) / mean(dplot$gene[dplot$ROI == "LC"])
mean(dplot$gene[dplot$ROI == "SNV"]) / mean(dplot$gene[dplot$ROI == "VTA"])

fact_lvls <- aggregate(gene ~ ROI, dplot, mean)$ROI[order(-aggregate(gene ~ ROI, dplot, mean)$gene)]
dplot$ROI <- factor(dplot$ROI, levels = fact_lvls)

dplot2 <- dplot[,c("ROI","gene")]
v1 <- violin_plot_function(dplot2,"ROI","gene","ROI",paste0(gene_name, " expression"), ROI_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="CTR ROI's (TH)")

arrange <- ggarrange(plotlist=list(v1), nrow=2, ncol=2)
ggsave("Violin_CTR.TH.ROI.SNCA.png", arrange)

# ii) high/ low SNCA expression cell-types
dplot2 <- dplot[dplot$ROI == "SNV",]
dplot2$high_SNCA <- dplot2$gene > 250

# get cell-type proportions and subset with 
myRCTD <- readRDS("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/RCTD_results/myRCTDde_gx.rds")
norm_weights <- normalize_weights(myRCTD@results$weights)
dplot3 <- merge(dplot2,norm_weights, by = "row.names")

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

# boxplot function
boxplot_2grp_function <- function(data_table, x_variable, y_variable, x_lab, y_lab, grp_variable,colour_palette) {

  # format data
  data_table$g <- data_table[, grp_variable]
  data_table$g <- as.factor(data_table$g)
  data_table$group <- factor(rep(c("grp1", "grp2"), nrow(data_table)/2))
  
  # order factor levels
  # perform t-test
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  data_table$g <- data_table[, grp_variable]
  
  stat.test <- data_table %>%
    group_by(x) %>%
    t_test(y ~ g) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  
  stat.test <- stat.test[order(unlist(stat.test[,"statistic"])),]
  fact_lvls <- stat.test$x[order(unlist(stat.test[,"statistic"]))]
  print(fact_lvls)
  data_table[, x_variable] <- factor(data_table[, x_variable] , levels = fact_lvls)
  
  # plot
  bxp <- ggboxplot(
    data_table, x = x_variable, y = y_variable,
    color = grp_variable, palette = colour_palette
  ) + ylab(y_lab) + xlab(x_lab)
   
  # add stat tests
  bxp + geom_pwc(
    aes(group = g), tip.length = 0,
    method = "t_test", label = "p.adj.signif", p.adjust.method = "BH"
  )
}

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




## scRNAseq ##
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


## fold changes



# ii) Violin plots between CTR & PD for each ROI


roi <- unique(dplot$ROI)
res <- list()
for (i in 1:length(roi)){
  print(roi[i])
  dplot2 <- dplot[dplot$segment == "TH" & dplot$ROI == roi[i],c("Diagnosis","gene")]

  res[[i]] <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
    theme(plot.title = element_text(hjust = 1)) +
    labs(title=roi[i])
}

arrange <- ggarrange(plotlist=list(v1,v2,v3,v4,v5,v6,v7,v8), nrow=2, ncol=2, 
                     labels=c("A9:Full ROI","A9:TH","SNV:Full ROI","SNV:TH"),
                     "SND:Full ROI","SND:TH","SNM:Full ROI","SNM:TH")

arrange <- ggarrange(plotlist=res, nrow=4, ncol=2)
ggsave("Violin_CTR.TH.ROI.SNCA.png", arrange)



## single-cell data ##


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


#### 
load(rdata)
meta_dat <- as.data.frame(gxdat_s@meta.data)
exp_dat <- as.matrix(gxdat_s@assays$RNA$counts)
tmp <- exp_dat[row.names(exp_dat) %in% NM_synthesis_genes,meta_dat$segment == "TH" & meta_dat$Diagnosis != "CTR"]
tmp2 <- meta_dat[meta_dat$segment == "TH" & meta_dat$Diagnosis != "CTR",]

cell_type <- as.data.frame(tmp2$ROI)
data_table <- cbind(cell_type,t(tmp))

# create plot data
dplot <- as.data.frame(data_table) %>%
  group_by(`tmp2$ROI`) %>% 
  summarise_all("median")
dplot1 <- as.matrix(t(dplot[,2:ncol(dplot)]))

# column annotations
Cell_Broad = dplot$`tmp2$ROI`
anno_df = data.frame(
  Region = dplot$`tmp2$ROI`
)

ha = HeatmapAnnotation(df = anno_df)

# # draw violin
# pdf("heatmap_snrnaseq_channelmembers_median_annot.pdf", width=15, height=15)
hm <- Heatmap(dplot1, cluster_columns = FALSE,
              top_annotation = ha,
              column_split=as.factor(Cell_Broad)
)
draw(hm,
     column_title = "NM synthesis genes in ST ILBD/ePD/lPD",
     column_title_gp=grid::gpar(fontsize=16))


log10(dplot1[row.names(dplot1) == "TH",]) / log10(dplot1[row.names(dplot1) == "SLC18A2",])
                                                                

log10(colSums(dplot1[row.names(dplot1) != "SLC18A2",])) / log10(dplot1[row.names(dplot1) == "SLC18A2",])

                                                                



### ST data - cell-type v SNCA levels
dat <- as.matrix(gxdat_s@assays$RNA$data)
gene_name <- "SNCA"
gene <- dat[gene_name,]
norm_weights <- normalize_weights(myRCTD@results$weights)
colnames(norm_weights) <- gsub("\\-","\\_",colnames(norm_weights))

dplot <- cbind(meta_dat,gene,norm_weights)

roi <- unique(dplot$ROI)
res <- list()
for (i in 1:length(roi)){
  print(roi[i])
  d1 <- dplot[dplot$ROI == roi[i] & dplot$segment == "TH",]
  model <- lm(formula(paste("gene ~ ",paste0(colnames(norm_weights),collapse = " + "))) , data = d1)
  res[[i]] <- summary(model)
}


############################################################################################
###### Part 2. Printing gene expression of interest
############################################################################################
## GeoMx data ##
load(rdata)
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
###### Part 3. Violin plots gene expression of interest
############################################################################################
# source plot functions
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/plotting_functions.R")

# geomx data
load(rdata)
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


