# Differential Expression Analysis
library(edgeR)
library(fgsea)
library(reactome.db)
library(metafor)
library(openxlsx)
library(viridis)
library(ggpubr)
library(dplyr)

############################################################################################
#### Inputs
############################################################################################
run_name = "geomx_oct2023"
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/"
rdata <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata"
contrast_path <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/contrasts_matrix.xlsx" # file with all regression models
NM_data <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/nm_data_090224.Rdata"

# source plotting functions
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/plotting_functions.R")

### Prepare data matrix, meta data and contrast matrices
# NOTE: The following scripts is designed to test list of contrasts (in contrasts_matrix.xlsx) by LIMMA Voom DEG analysis

############################################################################################
###### Part 1: LIMMA Voom
############################################################################################
# setwd
setwd(analysis_dir)
# Load normalized and batch corrected rdata
load(rdata)
# extract expression values
exp_dat <- as.matrix(gxdat_s@assays$RNA$counts)

# extract meta data
meta_dat <- as.data.frame(gxdat_s@meta.data)
table(row.names(meta_dat) == colnames(exp_dat))
meta_dat$scan_id <- row.names(meta_dat)

# format numerical variables
meta_dat$DV200 <- as.numeric(meta_dat$DV200)
meta_dat$Age <- as.numeric(meta_dat$Age)
meta_dat$n_per.µm2.iNM <- as.numeric(meta_dat$n_per.µm2.iNM)
meta_dat$n_per.µm2.eNM <- as.numeric(meta_dat$n_per.µm2.eNM)
meta_dat$Diagnosis_stage <- as.character(meta_dat$Diagnosis)
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "CTR"] <- 0
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "ILBD"] <- 1
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "ePD"] <- 2
meta_dat$Diagnosis_stage[meta_dat$Diagnosis_stage == "lPD"] <- 3
meta_dat$Diagnosis_stage <- as.numeric(meta_dat$Diagnosis_stage)
meta_dat$id <-  paste(meta_dat$Brainbank_ID,sep=".",meta_dat$ROI)

# load NM data 
load(NM_data)

# 1) summarise iNM Area per per SNV and LC per subject
data_table <- df_agg_iNM[df_agg_iNM$intra.extra == "iNM" ,c("Brainbank_ID","ROI","Area..µm..")]
df <- data_table %>%
  group_by(Brainbank_ID, ROI) %>%
  summarize_at(vars(-group_cols()), mean, na.rm = TRUE)
df$id <- paste(df$Brainbank_ID,sep=".",df$ROI)
# merge
meta_dat <- merge(meta_dat,df,by=c("id"), all = T)
meta_dat <- meta_dat[!is.na(meta_dat$`scan name`),]
colnames(meta_dat) <- make.names(colnames(meta_dat))
meta_dat$Brainbank_ID <- meta_dat$Brainbank_ID.x
meta_dat$ROI <- meta_dat$ROI.x
row.names(meta_dat) <- meta_dat$scan_id
# assign ROIs to high/ low iNM areas based on quantile's
quantiles <- quantile(meta_dat$Area..µm.., probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
meta_dat$iNM_quantile <- cut(meta_dat$Area..µm.., breaks = quantiles, include.lowest = TRUE, 
                       labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
meta_dat$iNM_quantile <- factor(meta_dat$iNM_quantile, levels = c("3rd_Q","1st_Q", "2nd_Q", "4th_Q") )



# 2) calculate group-wise % of toxic, protective etc iNM granules
data_table <- df_agg_iNM[!df_agg_iNM$iNM_cluster %in% "young CTR",]
colnames(data_table) <- make.names(colnames(data_table))
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID, ROI, Brainregion,Age,Sex,PMD,iNM_cluster,ROI.Area..µm.., sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$case.region <- paste0(dplot$Brainbank_ID,sep=".",dplot$ROI)

# select each clusetr
d1 <- dplot[dplot$iNM_cluster == "protective",]
d2 <- dplot[dplot$iNM_cluster == "novel",]
d3 <- dplot[dplot$iNM_cluster == "toxic",]
d4 <- dplot[dplot$iNM_cluster == "neutral",]

tmp <- merge(d1,d2[,c("case.region","n")], by="case.region",all=TRUE)
colnames(tmp)[colnames(tmp) %in% c("n.x","n.y")] <- c("n.1","n.2")
tmp <- merge(tmp,d3[,c("case.region","n")], by="case.region",all=TRUE)
colnames(tmp)[colnames(tmp) %in% c("n")] <- c("n.3")
tmp <- merge(tmp,d4[,c("case.region","n")], by="case.region",all=TRUE)
colnames(tmp)[colnames(tmp) %in% c("n")] <- c("n.4")

tmp[is.na(tmp)] <- 0
tmp <- tmp[complete.cases(tmp),]

tmp$protective_prc.total <- tmp$n.1 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
tmp$novel_prc.total <- tmp$n.2 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
tmp$toxic_prc.total <- tmp$n.3 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
tmp$neutral_prc.total <- tmp$n.4 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100

# assign to quantile's to SNV and LC separately 
tmp2 <- tmp[tmp$ROI == "SNV",]

quantiles <- quantile(tmp2$protective_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$protective_quantile<- cut(tmp2$protective_prc.total, breaks = quantiles, include.lowest = TRUE, 
                               labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$protective_quantile <- as.character(tmp2$protective_quantile)

quantiles <- quantile(tmp2$toxic_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$toxic_quantile<- cut(tmp2$toxic_prc.total, breaks = quantiles, include.lowest = TRUE, 
                          labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$toxic_quantile <- as.character(tmp2$toxic_quantile)

quantiles <- quantile(tmp2$novel_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$novel_quantile<- cut(tmp2$novel_prc.total, breaks = quantiles, include.lowest = TRUE, 
                          labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$novel_quantile <- as.character(tmp2$novel_quantile)

quantiles <- quantile(tmp2$neutral_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$neutral_quantile<- cut(tmp2$neutral_prc.total, breaks = quantiles, include.lowest = TRUE, 
                            labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$neutral_quantile <- as.character(tmp2$neutral_quantile)


tmp2 <- tmp[tmp$ROI == "LC",]

quantiles <- quantile(tmp2$protective_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$protective_quantile<- cut(tmp2$protective_prc.total, breaks = quantiles, include.lowest = TRUE, 
                               labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$protective_quantile <- as.character(tmp2$protective_quantile)

quantiles <- quantile(tmp2$toxic_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$toxic_quantile<- cut(tmp2$toxic_prc.total, breaks = quantiles, include.lowest = TRUE, 
                          labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$toxic_quantile <- as.character(tmp2$toxic_quantile)

quantiles <- quantile(tmp2$novel_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$novel_quantile<- cut(tmp2$novel_prc.total, breaks = quantiles, include.lowest = TRUE, 
                          labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$novel_quantile <- as.character(tmp2$novel_quantile)

quantiles <- quantile(tmp2$neutral_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$neutral_quantile<- cut(tmp2$neutral_prc.total, breaks = quantiles, include.lowest = TRUE, 
                            labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$neutral_quantile <- as.character(tmp2$neutral_quantile)




# mapvalues
meta_dat$protective_quantile <- mapvalues(meta_dat$id, from = tmp$case.region, to = tmp$protective_quantile)
meta_dat$protective_quantile[!meta_dat$protective_quantile %in% c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q")] <- 0

meta_dat$toxic_quantile<- mapvalues(meta_dat$id, from = tmp$case.region, to = tmp$toxic_quantile)
meta_dat$toxic_quantile[!meta_dat$toxic_quantile %in% c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q")] <- 0

meta_dat$novel_quantile <- mapvalues(meta_dat$id, from = tmp$case.region, to = tmp$novel_quantile)
meta_dat$novel_quantile[!meta_dat$novel_quantile %in% c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q")] <- 0

meta_dat$neutral_quantile <- mapvalues(meta_dat$id, from = tmp$case.region, to = tmp$neutral_quantile)
meta_dat$neutral_quantile[!meta_dat$neutral_quantile %in% c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q")] <- 0

# format factors
factor_format <- c("Brainbank_ID","Sex","Diagnosis","Brainregion","ROI")
for (i in 1:length(factor_format)){
  meta_dat[,factor_format[i]] <- as.factor(meta_dat[,factor_format[i]])
}
#meta_dat$Brainregion <- factor(meta_dat$Brainregion, levels=c('RN','A6','A9','A10'))
meta_dat$Diagnosis <- factor(meta_dat$Diagnosis, levels=c('ePD','CTR','lPD','ILBD'))

# read in contrast matrix
cont_dat <- read_excel(contrast_path, sheet = "LIMMA_Voom")
cont_dat <- cont_dat[cont_dat$notes == "run",]

## Run Voom
res <- list() # list to collect results
for (val in 147:nrow(cont_dat)){
  # print model number
  print(cont_dat$model.number[val])
  
  # create targets df
  if (cont_dat$roi[val] != "NA"){
    print("segment & roi")
    targ <- meta_dat[meta_dat$segment %in% unlist(strsplit(cont_dat$segment[val],",")) & meta_dat$ROI %in% unlist(strsplit(cont_dat$roi[val],",")),]
  } else if (cont_dat$brainregion[val] != "NA"){
    print("segment & brainregion")
    targ <- meta_dat[meta_dat$segment %in% unlist(strsplit(cont_dat$segment[val],",")) & meta_dat$Brainregion %in% unlist(strsplit(cont_dat$brainregion[val],",")),]
  } else {
    print("segment")
    targ <- meta_dat[meta_dat$segment %in% cont_dat$segment[val],]
  }
  
  # create contrast arg list
  cont_list <- list()
  cont_vars <- factor_format
  for (z in 1:length(cont_vars)){
    contrasts(targ[,cont_vars[z]], contrasts = FALSE)
    cont_list[[z]] <- contrasts(targ[,cont_vars[z]], contrasts = FALSE)
  }
  names(cont_list) <- factor_format
  
  # create design matrix
  options(na.action='na.omit')
  design <- model.matrix(reformulate(cont_dat$model.matrix[val]),
                         data=targ, 
                         drop = FALSE,
                         contrasts.arg = cont_list)
  
  # make names
  colnames(design) <- make.names(colnames(design))
  design <- design[,!colnames(design) %in% c("DiagnosisILBD.n_per.µm2.iNM")]
  
  # create expression df & targ with design row.names
  targ <- targ[row.names(design),]
  y <- exp_dat[,row.names(targ)]
  
  # fit design
  v <- voom(y,design)
  vfit <- lmFit(v)
  
  # Perform LIMMA contrasts
  cont.matrix <- makeContrasts(A=cont_dat$contrast[val],levels=design)
  # cont.matrix <- makeContrasts(A="DiagnosisCTR - DiagnosisILBD",levels=design)
  # cont.matrix <- makeContrasts(A="ROIRN",levels=design)
  fit2 <- contrasts.fit(vfit, cont.matrix)
  vfit2 <- eBayes(fit2)
  options(digits=3)
  
  # Select significant DEGs and assign to list
  tmp <- topTable(vfit2,number = Inf,sort.by="P")
  tmp2 <- tmp[tmp$adj.P.Val < 0.05,]
  res[[val]] <- tmp[tmp$P.Value < 0.05,]

}



# DEG enrichment within pigment related gene-sets
yf_pigment_genes <- unlist(read.delim(file="/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min/pigmentation+YuHong.genes.txt", header = F)[,1])

index1 <- which(names(vfit2$t[,1]) %in% yf_pigment_genes)
geneSetTest <- cameraPR(vfit2$t[,1],index1)

# i. name each item by description
voom_res <- res
names_1 <- paste0(cont_dat$description, "(",cont_dat$segment,")")
names(voom_res) <- names_1
for (i in 1:length(voom_res)){
  voom_res[[i]]$Contrast <- names(voom_res)[i]
}

# ## SAVE results list
# save(voom_res,cont_dat,file=paste0(analysis_dir,"voom_deg.rds"))

# write gene names
res <- voom_res
#res <- lapply(voom_res, function(x) x[x$adj.P.Val < 0.05,])
res <- lapply(res,function(x){
  x$Gene <- row.names(x)
  return(x)
} )

# write to excell by concept
diag <- res[which(cont_dat$concept_bin == "Diagnosis")]
write.xlsx(diag, file = "Diagnosis_LIMMA_Voom.xlsx", row.names = FALSE)

diag <- res[which(cont_dat$concept_bin == "Brainregion")]
write.xlsx(diag, file = "Brainregion_LIMMA_Voom.xlsx", row.names = FALSE)

diag <- res[which(cont_dat$concept_bin == "Neuromelanin")]
write.xlsx(diag, file = "Neuromelanin_LIMMA_Voom.xlsx", row.names = FALSE)

############################################################################################
###### Part 2: Exploratory plots of LIMMA Voom results
############################################################################################
## Barplots of n sig genes
## Neuromelanin in CTRs
cont_dat$contrast_short <- paste(paste0(cont_dat$description," [",cont_dat$segment,"]"))
cont_dat$DEG_adj.p0.05 <- as.numeric(cont_dat$DEG_adj.p0.05)

dplot <- cont_dat[cont_dat$concept_bin == "Neuromelanin",]
dplot <- dplot[!grepl("PD",dplot$contrast_short),]
g1 <- ggplot(dplot, aes(x=reorder(contrast_short,-DEG_adj.p0.05),y=DEG_adj.p0.05)) +
  geom_bar(stat="identity", width=0.7)+  
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("Contrast") + ylab("Sig DEG's (n)") 

dplot <- cont_dat[cont_dat$concept_bin == "Neuromelanin",]
dplot <- dplot[grepl("PD",dplot$contrast_short),]
g2 <- ggplot(dplot, aes(x=reorder(contrast_short,-DEG_adj.p0.05),y=DEG_adj.p0.05)) +
  geom_bar(stat="identity", width=0.7)+  
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("Contrast") + ylab("Sig DEG's (n)") 


### Select DEGs of diagnosis
# format cont_dat
contrast_short <- paste(paste0(cont_dat$brainregion," [",cont_dat$segment,"]"))
contrast_short[cont_dat$brainregion == "NA"] <- paste(paste0(cont_dat$roi," [",cont_dat$segment,"]"))[cont_dat$brainregion == "NA"]
contrast_short[cont_dat$brainregion == "NA" & cont_dat$roi == "NA"] <- paste0("[",cont_dat$segment,"]")[cont_dat$brainregion == "NA" & cont_dat$roi == "NA"] 
cont_dat$contrast_short <- paste0(contrast_short," - ",cont_dat$contrast)

cont_dat$DEG_adj.p0.05 <- as.numeric(cont_dat$DEG_adj.p0.05)

# barplot for n sig genes x contrast

dplot <- cont_dat[cont_dat$concept_bin == "Diagnosis",]
g2 <- ggplot(dplot, aes(x=reorder(contrast_short,-DEG_adj.p0.05),y=DEG_adj.p0.05)) +
  geom_bar(stat="identity", width=0.7)+  
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("Contrast") + ylab("Sig DEG's (n)") 

# save plots
ggsave("barplot_deg_dx_sigN.png", 
       ggarrange(g2,nrow=2,ncol=1)
       , device = "png")


# plots to check results by group
# example plots
dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat$Diagnosis<- factor(meta_dat$Diagnosis, levels=c('CTR','ILBD','ePD','lPD'))

gene_name <- "RAD23B"
gene <- dat[gene_name,]
dplot <- cbind(meta_dat,gene)
dplot <- dplot[dplot$segment == "TH",c("Diagnosis","gene")]

v1 <- violin_plot_function(dplot,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col)

gene_name <- "ACP3"
gene <- dat[gene_name,]
dplot <- cbind(meta_dat,gene)
dplot <- dplot[dplot$segment == "TH",c("Diagnosis_stage","gene")]

v2 <- violin_plot_function(dplot,"Diagnosis_stage","gene","Diagnosis Stage",paste0(gene_name, " expression"), Diagnosis_col)

arrange <- ggarrange(plotlist=list(v1,v2), nrow=2, ncol=2, widths = c(2,2))
ggsave("limma_voom_examples_dx.png", arrange)


# plots for genes of interest
gene_name <- "SNCA"
gene <- dat[gene_name,]
dplot <- cbind(meta_dat,gene)
dplot <- dplot[dplot$segment == "Full ROI" ,c("Diagnosis_stage","gene")]

violin_plot_function(dplot,"Diagnosis_stage","gene","Diagnosis Stage",paste0(gene_name, " expression"), Diagnosis_col)

gene_name <- "RGN"
gene <- log10(dat[gene_name,])
dplot <- cbind(meta_dat,gene)
dplot <- dplot[dplot$segment == "Full ROI" ,c("ROI","gene")]
dplot$ROI <- factor(dplot$ROI, levels = levels(reorder(dplot$ROI, dplot$gene, mean)))

violin_plot_function_p(dplot,"ROI","gene","ROI",paste0(gene_name, " expression"), ROI_col)

gene_name <- "ZFP36"
gene <- log10(dat[gene_name,])
dplot <- cbind(meta_dat,gene)
dplot <- dplot[dplot$segment == "Full ROI"  ,]
ggplot(dplot, aes(x=n.iNM, y=gene,  color=Diagnosis)) +
  geom_point() +
  geom_smooth(method=lm, se =FALSE)

lapply(voom_res, function(x) x[row.names(x) == "RGN",])
lapply(cside_res, function(x) {
  x[row.names(x) == "RGN",]
  
  }
  )
