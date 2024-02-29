# Differential Expression Analysis
library(edgeR)
library(fgsea)
library(reactome.db)
library(metafor)
library(openxlsx)
library(viridis)
library(ggpubr)
library(dplyr)
library(plyr)

############################################################################################
#### Inputs
############################################################################################
run_name = "geomx_oct2023"
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/"
rdata <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata"
contrast_path <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/contrasts_matrix.xlsx" # file with all regression models
#NM_data <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/nm_data_090224.Rdata"
NM_data <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/nm_data_280224.Rdata"

# source plotting functions
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/plotting_functions.R")

### Prepare data matrix, meta data and contrast matrices
# NOTE: The following scripts is designed to test list of contrasts (in contrasts_matrix.xlsx) by LIMMA Voom DEG analysis

# colour palettes
Diagnosis_col = c("CTR"= "grey","ILBD" = "#00AFBB", "ePD" = "#E7B800","lPD" = "red")
age_group_col = magma(4)
ROI_col = c("SNV" = "purple",
            "SNM" = viridis(6)[2],
            "SND" = viridis(6)[3],
            "SNL" = viridis(6)[4],
            "VTA" = viridis(6)[5],
            "LC" = viridis(6)[6])

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
meta_dat$GenesDetected <- as.numeric(meta_dat$GenesDetected)

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



# 2) calculate group-wise % of toxic, nontoxic etc iNM granules
data_table <- df_agg_iNM[!df_agg_iNM$iNM_cluster %in% "young CTR",]
colnames(data_table) <- make.names(colnames(data_table))
dplot <- data_table %>% dplyr::count(intra.extra,Diagnosis,Brainbank_ID, ROI, Brainregion,Age,Sex,PMD,iNM_cluster,ROI.Area..µm.., sort = TRUE)
dplot <- dplot[!is.na(dplot$intra.extra),]
dplot$case.region <- paste0(dplot$Brainbank_ID,sep=".",dplot$ROI)

# select each clusetr
d1 <- dplot[dplot$iNM_cluster == "non-toxic",]
d2 <- dplot[dplot$iNM_cluster == "transition",]
d3 <- dplot[dplot$iNM_cluster == "toxic",]
d4 <- dplot[dplot$iNM_cluster == "novel",]

tmp <- merge(d1,d2[,c("case.region","n")], by="case.region",all=TRUE)
colnames(tmp)[colnames(tmp) %in% c("n.x","n.y")] <- c("n.1","n.2")
tmp <- merge(tmp,d3[,c("case.region","n")], by="case.region",all=TRUE)
colnames(tmp)[colnames(tmp) %in% c("n")] <- c("n.3")
tmp <- merge(tmp,d4[,c("case.region","n")], by="case.region",all=TRUE)
colnames(tmp)[colnames(tmp) %in% c("n")] <- c("n.4")

tmp[is.na(tmp)] <- 0
tmp <- tmp[complete.cases(tmp),]

tmp$nontoxic_prc.total <- tmp$n.1 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
tmp$transition_prc.total <- tmp$n.2 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
tmp$toxic_prc.total <- tmp$n.3 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100
tmp$novel_prc.total <- tmp$n.4 / (tmp$n.1 + tmp$n.2 + tmp$n.3 + tmp$n.4) * 100

# assign to quantile's to SNV and LC separately 
tmp2 <- tmp[tmp$ROI == "SNV",]

quantiles <- quantile(tmp2$nontoxic_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$nontoxic_quantile<- cut(tmp2$nontoxic_prc.total, breaks = quantiles, include.lowest = TRUE, 
                               labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$nontoxic_quantile <- as.character(tmp2$nontoxic_quantile)

quantiles <- quantile(tmp2$toxic_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$toxic_quantile<- cut(tmp2$toxic_prc.total, breaks = quantiles, include.lowest = TRUE, 
                          labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$toxic_quantile <- as.character(tmp2$toxic_quantile)

quantiles <- quantile(tmp2$novel_prc.total, probs = c(0, 0.5, 0.75, 1),na.rm =T)
tmp2$novel_quantile<- cut(tmp2$novel_prc.total, breaks = quantiles, include.lowest = TRUE, 
                          labels = c("1st.2nd_Q", "3rd_Q", "4th_Q")) # note labelled as 1st_Q for regression purposed
tmp2$novel_quantile <- as.character(tmp2$novel_quantile)

quantiles <- quantile(tmp2$transition_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$transition_quantile<- cut(tmp2$transition_prc.total, breaks = quantiles, include.lowest = TRUE, 
                            labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$transition_quantile <- as.character(tmp2$transition_quantile)

tmp_snv <- tmp2

tmp2 <- tmp[tmp$ROI == "LC",]

quantiles <- quantile(tmp2$nontoxic_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$nontoxic_quantile<- cut(tmp2$nontoxic_prc.total, breaks = quantiles, include.lowest = TRUE, 
                             labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$nontoxic_quantile <- as.character(tmp2$nontoxic_quantile)

quantiles <- quantile(tmp2$toxic_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$toxic_quantile<- cut(tmp2$toxic_prc.total, breaks = quantiles, include.lowest = TRUE, 
                          labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$toxic_quantile <- as.character(tmp2$toxic_quantile)

quantiles <- quantile(tmp2$transition_prc.total, probs = c(0, 0.25, 0.5, 0.75, 1),na.rm =T)
tmp2$transition_quantile<- cut(tmp2$transition_prc.total, breaks = quantiles, include.lowest = TRUE, 
                               labels = c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q"))
tmp2$transition_quantile <- as.character(tmp2$transition_quantile)
tmp2$novel_quantile <- NA

tmp <- rbind(tmp_snv,tmp2)

# mapvalues
meta_dat$nontoxic_quantile <- mapvalues(meta_dat$id, from = tmp$case.region, to = tmp$nontoxic_quantile)
meta_dat$nontoxic_quantile[!meta_dat$nontoxic_quantile %in% c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q")] <- 0

meta_dat$toxic_quantile<- mapvalues(meta_dat$id, from = tmp$case.region, to = tmp$toxic_quantile)
meta_dat$toxic_quantile[!meta_dat$toxic_quantile %in% c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q")] <- 0

meta_dat$novel_quantile <- mapvalues(meta_dat$id, from = tmp$case.region, to = tmp$novel_quantile)
meta_dat$novel_quantile[!meta_dat$novel_quantile %in% c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q")] <- 0

meta_dat$transition_quantile <- mapvalues(meta_dat$id, from = tmp$case.region, to = tmp$transition_quantile)
meta_dat$transition_quantile[!meta_dat$transition_quantile %in% c("1st_Q", "2nd_Q", "3rd_Q", "4th_Q")] <- 0

# plot the overlap in ROI's
library(VennDiagram)

venn.diagram(
  x = list(row.names(meta_dat)[meta_dat$nontoxic_quantile == "1st_Q"],
           row.names(meta_dat)[meta_dat$toxic_quantile == "1st_Q"], 
           row.names(meta_dat)[meta_dat$novel_quantile == "1st_Q"], 
           row.names(meta_dat)[meta_dat$transition_quantile == "1st_Q"]),
  category.names = c("nontoxic" , "Toxic" , "Novel","transition"),
  filename = 'ROI_iNMQ1.2_venn_diagramm.png',
  output=TRUE,main = "Q1.2"
)

venn.diagram(
  x = list(row.names(meta_dat)[meta_dat$nontoxic_quantile == "4th_Q"],
           row.names(meta_dat)[meta_dat$toxic_quantile == "4th_Q"], 
           row.names(meta_dat)[meta_dat$novel_quantile == "4th_Q"], 
           row.names(meta_dat)[meta_dat$transition_quantile == "4th_Q"]),
  category.names = c("nontoxic" , "Toxic" , "Novel","transition"),
  filename = 'ROI_iNMQQ4_venn_diagramm.png',
  output=TRUE,main = "Q4"
)

venn.diagram(
  x = list(row.names(meta_dat)[meta_dat$novel_quantile == "1st_Q"], 
           row.names(meta_dat)[meta_dat$transition_quantile == "4th_Q"]),
  category.names = c("Novel","transition"),
  filename = 'ROI_iNM_novel1.transition4_venn_diagramm.png',
  output=TRUE,main = "novel1 transition4"
)

venn.diagram(
  x = list(row.names(meta_dat)[meta_dat$novel_quantile == "4th_Q"], 
           row.names(meta_dat)[meta_dat$transition_quantile == "1st_Q"]),
  category.names = c("Novel","transition"),
  filename = 'ROI_iNM_novel4.transition1_venn_diagramm.png',
  output=TRUE,main = "novel4 transition1"
)


# format factors
factor_format <- c("Brainbank_ID","Sex","Diagnosis","Brainregion","ROI")
for (i in 1:length(factor_format)){
  meta_dat[,factor_format[i]] <- as.factor(meta_dat[,factor_format[i]])
}
#meta_dat$Brainregion <- factor(meta_dat$Brainregion, levels=c('RN','A6','A9','A10'))
meta_dat$Diagnosis <- factor(meta_dat$Diagnosis, levels=c('ePD','CTR','lPD','ILBD'))

# read in contrast matrix
library(writexl)
cont_dat <- read_xlsx(contrast_path, sheet = "LIMMA_Voom")
cont_dat <- cont_dat[cont_dat$notes == "run",]

## Run Voom
res <- list() # list to collect results
sample_names_list <- list()
for (val in 146:nrow(cont_dat)){
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
  
  # store names of contrasts
  sample_names_list[[val]] <- row.names(fit2$design)
  
  # Select significant DEGs and assign to list
  tmp <- topTable(vfit2,number = Inf,sort.by="P")
  tmp2 <- tmp[tmp$adj.P.Val < 0.05,]
  print(nrow(tmp2))
  res[[val]] <- tmp

}


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
# diag <- res[which(cont_dat$concept_bin == "Diagnosis")]
# write.xlsx(diag, file = "Diagnosis_LIMMA_Voom.xlsx", row.names = FALSE)
# 
# diag <- res[which(cont_dat$concept_bin == "Brainregion")]
# write.xlsx(diag, file = "Brainregion_LIMMA_Voom.xlsx", row.names = FALSE)
# 
# diag <- res[which(cont_dat$concept_bin == "Neuromelanin")]
# write.xlsx(diag, file = "Neuromelanin_LIMMA_Voom.xlsx", row.names = FALSE)

diag <- res[146:length(res)]
write.xlsx(diag, file = "Neuromelanin_classified_LIMMA_Voom_v2.xlsx", row.names = FALSE)


############################################################################################
###### Part 2: Gene ontology
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

# create list 
gene_interest_list <- list(skin_melanin_enzymes,alt_PD,CA_precursor_NM_stock_trans_metab_A9.A10,CA_precursor_NM_stock_trans_metab_A6,
                           CA_functional_DA,CA_functional_NE,Stress_granules.free_radical_scavenging,Lysosome_pathways,upstream_euNM,
                           genes_link_skinpig.PD,pigmentation_network,other)

names(gene_interest_list) <- c("skin_melanin_enzymes","alt_PD","CA_precursor_NM_stock_trans_metab_A9.A10","CA_precursor_NM_stock_trans_metab_A6",
                              "CA_functional_DA","CA_functional_NE","Stress_granules.free_radical_scavenging","Lysosome_pathways","upstream_euNM",
                              "genes_link_skinpig.PD","pigmentation_network","other")

##  Identify ontologies for each contrast
cont_list <- res[146:length(res)]
# db
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
# run
up_ont_list <- list()
down_ont_list <- list()
for (i in 1:length(cont_list)){
  print(names(cont_list[i]))
  exp_cont <- as.data.frame(cont_list[[i]])
  
  # upreg list
  upreg <- exp_cont$Gene[exp_cont$adj.P.Val < 0.05 & exp_cont$logFC > 0]
  print(length(upreg))
  
  # downreg list
  downreg <- exp_cont$Gene[exp_cont$adj.P.Val < 0.05 & exp_cont$logFC < 0]
  print(length(downreg))

  # upreg top 3 passing adj p-val < 0.05
  enriched <- enrichr(upreg, dbs)
  enriched <- lapply(enriched, function(x) {
    x <- x[x$Adjusted.P.value < 0.05,]
  })
  filtered_list <- Filter(function(x) nrow(x) > 0, enriched)
  up_ont_list[[i]] <- bind_rows(filtered_list, .id = "List_Name")
  
  # downreg top 3 passing adj p-val < 0.05
  enriched <- enrichr(downreg, dbs)
  enriched <- lapply(enriched, function(x) {
    x <- x[x$Adjusted.P.value < 0.05,]
  })
  
  filtered_list <- Filter(function(x) nrow(x) > 0, enriched)
  down_ont_list[[i]] <- bind_rows(filtered_list, .id = "List_Name")
  
}



############################################################################################
###### Part 3: Heatmaps
############################################################################################
cont_val = as.numeric(which(lapply(voom_res,function(x) any(x$adj.P.Val < 0.05)) == "TRUE")) # contrasts with sig genes

for (rel_val in 1:length(cont_val)){
  print(rel_val)
  # get dataset from contrast
  val=cont_val[rel_val] # contrast numbers
  
  set.seed(123)
  # print model number
  print(cont_dat$model.number[val])
  cont_dat$description[val]
  
  # create targets df
  targ <- meta_dat[sample_names_list[[val]],]
  targ <- targ[targ[,unlist(strsplit(cont_dat$description[val],","))[1]] %in% c("1st_Q","1st.2nd_Q","4th_Q"),]
  
  # get sig genes
  genes <- voom_res[[val]]
  genes <- row.names(genes)[genes$adj.P.Val < 0.05]
  
  # create expression df & targ with design row.names
  y <- exp_dat[genes,row.names(targ)]
  
  # retrieve GO row annotations
  tmp <- up_ont_list[[rel_val]]
  tmp <- tmp[tmp$Adjusted.P.value < 0.05,]
  
  go_list <- list()
  tmp2 <- tmp[tmp$List_Name == "GO_Biological_Process_2015",][1:3,]
  if(length(tmp2) > 0){
    list_bool <- lapply(strsplit(tmp2$Genes,";"), function(x) genes %in% x)
    tmp3 <- t(do.call(rbind,list_bool))
    colnames(tmp3) <- tmp2$Term
    row.names(tmp3) <- genes
    go_list[[1]] <- tmp3
  }
  
  tmp2 <- tmp[tmp$List_Name == "GO_Cellular_Component_2015",][1:3,]
  if(length(tmp2) > 0){
    list_bool <- lapply(strsplit(tmp2$Genes,";"), function(x) genes %in% x)
    tmp3 <- t(do.call(rbind,list_bool))
    colnames(tmp3) <- tmp2$Term
    row.names(tmp3) <- genes
    go_list[[2]] <- tmp3
  }
  
  tmp2 <- tmp[tmp$List_Name == "GO_Molecular_Function_2015",][1:3,]
  if(length(tmp2) > 0){
    list_bool <- lapply(strsplit(tmp2$Genes,";"), function(x) genes %in% x)
    tmp3 <- t(do.call(rbind,list_bool))
    colnames(tmp3) <- tmp2$Term
    row.names(tmp3) <- genes
    go_list[[3]] <- tmp3
  }
  
  go_annot <- do.call(cbind,go_list)
  
  # compare to 
  list_bool <- lapply(gene_interest_list, function(x) genes %in% x)
  tmp3 <- t(do.call(rbind,list_bool))
  colnames(tmp3) <- names(gene_interest_list)
  row.names(tmp3) <- genes
  
  annotation_matrix <- cbind(go_annot,tmp3)
  annotation_matrix <- apply(annotation_matrix,2,as.factor)
  
  # # restrict to any annots
  # any_annot <- apply(annotation_matrix,1,any)
  # annotation_matrix <- annotation_matrix[any_annot,]
  # y <- y[any_annot,]
  
  # construct row annotations 
  row_ha <- rowAnnotation("Annotation" = annotation_matrix, col = list("TRUE" = "grey", "FALSE" = "black"), border = TRUE,gp = gpar(col = "black"))
  
  # construct column annotations
  col_ha = HeatmapAnnotation(
    Dx = targ$Diagnosis,
    ROI = targ$ROI,
    Q = targ[,unlist(strsplit(cont_dat$description[val],","))[1]],
    col = list(Dx = Diagnosis_col, ROI = ROI_col , Q = c("1st_Q" = "lightblue","1st.2nd_Q" = "lightblue","4th_Q" = "darkblue") )
  )
  
  # plot heatmap
  ht <- Heatmap(log10(y), name = "Log10(Counts)", 
                right_annotation = row_ha,
                top_annotation = col_ha, 
                show_column_names = FALSE,
                cluster_columns = TRUE,
                row_names_gp = gpar(fontsize = 3),
                #row_split = targets$Group,
                column_title = cont_dat$description[val],
                cluster_row_slices = FALSE,
  
                #row_title_gp = gpar(col = c("Dopamine" = "purple","Ca2+ & Transmission" = "orange"), cex = 0.7),
                 width = ncol(y)*unit(1, "mm"), 
                 height = nrow(y)*unit(1, "mm"))
  
  pdf(paste0("heatmap_",make.names(cont_dat$description[val]),".",val,".pdf"),width = 15,height =100)
  draw(ht)
  dev.off()
}


# assign plot data
name_plot = "hits_aavgq_genes.interest"
#name_plot = "hits_aavgq.SN.TH_n59"

# plot
pdf(paste0("heatmap_",name_plot,".pdf"))
draw(ht, annotation_legend_list = list(lgd_pvalue, lgd_sig))
dev.off()













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


############################################################################################
###### Part 3: Exploratory plots of LIMMA Voom results
############################################################################################
# read in 


# DEG enrichment within pigment related gene-sets
yf_pigment_genes <- unlist(read.delim(file="/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis_min/pigmentation+YuHong.genes.txt", header = F)[,1])

index1 <- which(names(vfit2$t[,1]) %in% yf_pigment_genes)
geneSetTest <- cameraPR(vfit2$t[,1],index1)



