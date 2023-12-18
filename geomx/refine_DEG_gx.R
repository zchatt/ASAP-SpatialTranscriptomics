# integrating DEG and ctDEG adn prioritising DEGs for validation

library(ggVennDiagram)
library(STRINGdb)
library(rbioapi)
library(enrichR)

###############
### Input
###############
# analysis dir
working_dir = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
setwd(working_dir)

# ST data
rdata <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata"

# source plotting functions and poublic degs
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/plotting_functions.R")
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/deg_markers/format_deg_markers.R")

# PD genetics
g4pd_cnv_path <- "/Users/zacc/USyd/spatial_transcriptomics/data/ref_genesets/t_cnv20200820.txt"
g4pd_cv_path <- "/Users/zacc/USyd/spatial_transcriptomics/data/ref_genesets/t_common_variant20200820.txt"
g4pd_rg_path <- "/Users/zacc/USyd/spatial_transcriptomics/data/ref_genesets/t_rare_gene20200820.txt"
g4pd_rv_path <- "/Users/zacc/USyd/spatial_transcriptomics/data/ref_genesets/t_rare_variant20200820.txt"
g4pd_deg_path <- "/Users/zacc/USyd/spatial_transcriptomics/data/ref_genesets/t_gene_expression20200820.txt"
contrast_path <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/contrasts_matrix.xlsx"
genetic_list <- list()
genetic_list[[1]] <- read.delim(g4pd_rg_path, sep='\t', row.names=NULL, header=T)
genetic_list[[2]] <- read.delim(g4pd_cnv_path, sep='\t', row.names=NULL, header=T)
genetic_list[[3]] <- read.delim(g4pd_cv_path, sep='\t', row.names=NULL, header=T)
genetic_list[[4]] <- read.delim(g4pd_rv_path, sep='\t', row.names=NULL, header=T)
g4pd_deg <- read.delim(g4pd_deg_path)
colnames(g4pd_deg) <- paste0("Published_PD_DEG_evidence_",colnames(g4pd_deg))

# DEG and ctDEG results
load("cside_deg.rds") # cside results
load("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/voom_deg.rds") # limma voom results

# contrast matrices
cont_dat <- read_excel(contrast_path, sheet = "LIMMA_Voom")
contl <- cont_dat[cont_dat$notes == "run",]
cont_dat <- read_excel(contrast_path, sheet = "C-SIDE")
contc <- cont_dat[cont_dat$notes == "run",]

#################
### Functions ###
#################
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

#########################################################
## 1. Identify positive and negative regulators of NM ###
#########################################################
ctr.epd.lpd.full <- as.data.frame(voom_res[[132]])
ctr.epd.lpd.th <- as.data.frame(voom_res[[133]])
ctr.ilbd.full <- as.data.frame(voom_res[[134]])
ctr.ilbd.th <- as.data.frame(voom_res[[135]])
ctr.vtaelse.full <- as.data.frame(voom_res[[13]])
ctr.vtaelse.th <- as.data.frame(voom_res[[19]])

### FULL ROI ###
# bind dfs get pos logFC
datp <- rbind(ctr.epd.lpd.full)
# bind dfs to get neg logFC
datn<- rbind(ctr.vtaelse.full,ctr.ilbd.full)

# positive regulators NM
pos_reg <- rbind(datp[datp$logFC > 0,], datn[datn$logFC < 0,])
pos_reg <- pos_reg[pos_reg$adj.P.Val < 0.05,]
# supported by multiple lines of evidence
tmp <- table(pos_reg$Gene)
pos_reg2 <- pos_reg[pos_reg$Gene %in% names(tmp)[tmp > 1], ] 

#plot
tab1 <- as.data.frame.matrix(table(pos_reg2$Gene,pos_reg2$Contrast))
x = list(
  PD.ILBD = row.names(tab1)[tab1$`CTRvsePDlPD,alltissue(Full ROI)` == 1 & tab1$`CTRvsILBD,alltissue(Full ROI)` == 1],
  PD.VTA = row.names(tab1)[tab1$`CTRvsePDlPD,alltissue(Full ROI)`== 1  & tab1$`VTAvsElse(Full ROI)` == 1 ],
  ILBD.VTA = row.names(tab1)[tab1$`CTRvsILBD,alltissue(Full ROI)` == 1  & tab1$`VTAvsElse(Full ROI)` == 1 ]
  )
v1 <- ggVennDiagram(x, label_alpha = 0) + 
  scale_fill_distiller(palette = "Reds", direction = 1) +
  theme(legend.position = "none")

# negative regulators NM
neg_reg <- rbind(datp[datp$logFC < 0,],datn[datn$logFC > 0,])
neg_reg <- neg_reg[neg_reg$adj.P.Val < 0.05,]
# supported by multiple lines of evidence
tmp <- table(neg_reg$Gene)
neg_reg2 <- neg_reg[neg_reg$Gene %in% names(tmp)[tmp > 1], ] 

# plot
tab1 <- as.data.frame.matrix(table(neg_reg2$Gene,neg_reg2$Contrast))
x = list(
  PD.ILBD = row.names(tab1)[tab1$`CTRvsePDlPD,alltissue(Full ROI)` == 1 & tab1$`CTRvsILBD,alltissue(Full ROI)` == 1],
  PD.VTA = row.names(tab1)[tab1$`CTRvsePDlPD,alltissue(Full ROI)`== 1  & tab1$`VTAvsElse(Full ROI)` == 1 ],
  ILBD.VTA = row.names(tab1)[tab1$`CTRvsILBD,alltissue(Full ROI)` == 1  & tab1$`VTAvsElse(Full ROI)` == 1 ]
)
v2 <- ggVennDiagram(x, label_alpha = 0) + 
  scale_fill_distiller(palette = "Reds", direction = 1) +
  theme(legend.position = "none")

# cannot be in both - None are in both
neg_reg_full <- neg_reg2[!neg_reg2$Gene %in% pos_reg2$Gene,]
pos_reg_full <- pos_reg2[!pos_reg2$Gene %in% neg_reg2$Gene,]

# evaluating known NM synthesis members
NM_synthesis_genes <- c("TH","TYR","DDC","COMT","MAOA","MAOB","AKR1B1","SLC18A2","RGN",row.names(sc_obj)[grep("ALDH",row.names(sc_obj))])
nm_pos <- pos_reg_full[pos_reg_full$Gene %in% NM_synthesis_genes,]
nm_neg <- neg_reg_full[neg_reg_full$Gene %in% NM_synthesis_genes,]


### TH masked ROI ###
# bind dfs get pos logFC
datp <- rbind(ctr.epd.lpd.th)
# bind dfs to get neg logFC
datn<- rbind(ctr.vtaelse.th,ctr.ilbd.th)

# positive regulators NM
pos_reg <- rbind(datp[datp$logFC > 0,], datn[datn$logFC < 0,])
pos_reg <- pos_reg[pos_reg$adj.P.Val < 0.05,]
# supported by multiple lines of evidence
tmp <- table(pos_reg$Gene)
pos_reg2 <- pos_reg[pos_reg$Gene %in% names(tmp)[tmp > 1], ] 

#plot
tab1 <- as.data.frame.matrix(table(pos_reg2$Gene,pos_reg2$Contrast))
x = list(
  PD.ILBD = row.names(tab1)[tab1$`CTRvsePDlPD,alltissue(TH)` == 1 & tab1$`CTRvsILBD,alltissue(TH)` == 1],
  PD.VTA = row.names(tab1)[tab1$`CTRvsePDlPD,alltissue(TH)`== 1  & tab1$`VTAvsElse(TH)` == 1 ],
  ILBD.VTA = row.names(tab1)[tab1$`CTRvsILBD,alltissue(TH)` == 1  & tab1$`VTAvsElse(TH)` == 1 ]
)
v3 <- ggVennDiagram(x, label_alpha = 0) + 
  scale_fill_distiller(palette = "Reds", direction = 1)  +
  theme(legend.position = "none")


# negative regulators NM
neg_reg <- rbind(datp[datp$logFC < 0,],datn[datn$logFC > 0,])
neg_reg <- neg_reg[neg_reg$adj.P.Val < 0.05,]
# supported by multiple lines of evidence
tmp <- table(neg_reg$Gene)
neg_reg2 <- neg_reg[neg_reg$Gene %in% names(tmp)[tmp > 1], ] 

# plot
tab1 <- as.data.frame.matrix(table(neg_reg2$Gene,neg_reg2$Contrast))
x = list(
  PD.ILBD = row.names(tab1)[tab1$`CTRvsePDlPD,alltissue(TH)` == 1 & tab1$`CTRvsILBD,alltissue(TH)` == 1],
  PD.VTA = row.names(tab1)[tab1$`CTRvsePDlPD,alltissue(TH)`== 1  & tab1$`VTAvsElse(TH)` == 1 ],
  ILBD.VTA = row.names(tab1)[tab1$`CTRvsILBD,alltissue(TH)` == 1  & tab1$`VTAvsElse(TH)` == 1 ]
)
v4 <- ggVennDiagram(x, label_alpha = 0) + 
  scale_fill_distiller(palette = "Reds", direction = 1) +
  theme(legend.position = "none")

# cannot be in both - None are in both
neg_reg_th <- neg_reg2[!neg_reg2$Gene %in% pos_reg2$Gene,]
pos_reg_th <- pos_reg2[!pos_reg2$Gene %in% neg_reg2$Gene,]

# evaluating known NM synthesis members
NM_synthesis_genes <- c("TH","TYR","DDC","COMT","MAOA","MAOB","AKR1B1","SLC18A2","RGN",row.names(sc_obj)[grep("ALDH",row.names(sc_obj))])
nm_pos_th <- pos_reg_th[pos_reg_th$Gene %in% NM_synthesis_genes,]
nm_neg_th <- neg_reg_th[neg_reg_th$Gene %in% NM_synthesis_genes,]

# save plots
arrange <- ggarrange(plotlist=list(v1,v2,v3,v4), nrow=2, ncol=2, widths = c(2,2),labels=c("Full ROI Pos Reg","Full ROI Neg Reg","TH Pos Reg","TH Neg Reg"))
ggsave("venns_posneg_reg.png", arrange)


## format before writing
# list of pos and neg regulators
tmp_all <- list(neg_reg_th, neg_reg_full, pos_reg_th, pos_reg_full)

# calculate Stouffer p for 
for (z in 1:length(tmp_all)){
  tmp <- tmp_all[[z]]
  tmp$stouffer_p <- NA
  for (i in 1:nrow(tmp)){
    gene_1 <- tmp$Gene[i]
    p_val <- tmp$adj.P.Val[tmp$Gene == gene_1]
    stouffer_p <- combine_pvalues_stouffer(p_val)
    
    tmp$stouffer_p[tmp$Gene == gene_1] <- stouffer_p
  }
}

# Cross check published data
deg_df <- do.call(rbind, deg_list)# create dataframe from published genesets
tmp_all <- lapply(tmp_all, function(x) merge(x,deg_df,by="Gene", all.x = T))# merge

# Cross check published data
# create dataframe from published genesets
deg_df <- do.call(rbind, deg_list)
# merge
tmp_all <- lapply(tmp_all, function(x) merge(x,deg_df,by="Gene", all.x = T))

# Cross check genetic associations
gene_pd <- unlist(lapply(genetic_list, function(x) x$Gene))
# merge
tmp_all <- lapply(combined_df, function(x) {
  x$Genetic_Association <- x$Gene %in% gene_pd
  return(x)
})

# order
tmp_all <- lapply(tmp_all, function(x) {
  x <- x[order(x$stouffer_p),]
  return(x)
})

# write to excel
names(tmp_all) <- c("neg_reg_th", "neg_reg_full", "pos_reg_th", "pos_reg_full")
NM_reg_list <- tmp_all
write.xlsx(NM_reg_list, file = "Regulators_NM_DEG.xlsx", row.names = FALSE)

### Plots
# example plots
load(rdata)
dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat$Diagnosis<- factor(meta_dat$Diagnosis, levels=c('CTR','ILBD','ePD','lPD'))

# gene interest 
gene_name <- "CHCHD10"
gene <- dat[gene_name,]
dplot <- cbind(meta_dat,gene)

dplot2 <- dplot[dplot$segment == "Full ROI" & dplot$Brainregion == "A9",c("Diagnosis","gene")]
v1 <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="A9:Full ROI")

dplot2 <- dplot[dplot$segment == "TH" & dplot$Brainregion == "A9",c("Diagnosis","gene")]
v2 <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="A9:TH")

dplot2 <- dplot[dplot$segment == "Full ROI" & dplot$ROI == "SNV",c("Diagnosis","gene")]
v3 <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="SNV:Full ROI")

dplot2 <- dplot[dplot$segment == "TH" & dplot$ROI == "SNV",c("Diagnosis","gene")]
v4 <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="SNV:TH")

dplot2 <- dplot[dplot$segment == "Full ROI" & dplot$ROI == "SND",c("Diagnosis","gene")]
v5 <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="SND:Full ROI")

dplot2 <- dplot[dplot$segment == "TH" & dplot$ROI == "SND",c("Diagnosis","gene")]
v6 <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="SND:TH")

dplot2 <- dplot[dplot$segment == "Full ROI" & dplot$ROI == "SNM",c("Diagnosis","gene")]
v7 <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="SNM:Full ROI")

dplot2 <- dplot[dplot$segment == "TH" & dplot$ROI == "SNM",c("Diagnosis","gene")]
v8 <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="SNM:TH")


arrange <- ggarrange(plotlist=list(v1,v2,v3,v4,v5,v6,v7,v8), nrow=2, ncol=2, 
                     labels=c("A9:Full ROI","A9:TH","SNV:Full ROI","SNV:TH"),
                     "SND:Full ROI","SND:TH","SNM:Full ROI","SNM:TH")

arrange <- ggarrange(plotlist=list(v1,v2,v3,v4,v5,v6,v7,v8), nrow=4, ncol=2)
ggsave("Violin_CHCHD10_Dx.Region.ROI.png", arrange)


gene_name <- "KPNA3"
gene <- dat[gene_name,]
dplot <- cbind(meta_dat,gene)

dplot2 <- dplot[dplot$segment == "Full ROI" & dplot$ROI == "SND",c("Diagnosis","gene")]
v1 <- violin_plot_function(dplot2,"Diagnosis","gene","Diagnosis",paste0(gene_name, " expression"), Diagnosis_col) +
  theme(plot.title = element_text(hjust = 1)) +
  labs(title="SND:Full ROI")


# Get gene expression for Nakamura
gene_name <- c("CHCHD10","CHCHD2","SNCA","PDHA1")
gene <- dat[gene_name,]
dplot <- cbind(meta_dat,t(gene))
dplot <- as.data.frame(dplot)
dplot <- dplot[dplot$segment == "TH" & dplot$ROI %in% c("SNV","SND"),]
dplot <- dplot[,c("Brainbank_ID","Sex","Age","PMD hs","ROI","SNCA","PDHA1","CHCHD10","CHCHD2")]
write.xlsx(dplot, file = "targeted_genelist_geomx_011223.xlsx")


### Running STRING  ###
## Prioritizing NM genes by protein-protein interaction with Neuromelanin and Melanin Biosynthetic pathway members
NM_synthesis_genes_expressed <- NM_synthesis_genes[NM_synthesis_genes %in% row.names(dat)]
melanin_biosynthetic_process <- c("CITED1","TYR","PMEL","OCA2","MC1R","SLC24A5","DCT","MC1R","SLC45A2","PKS1") #GO:0042438
melanin_biosynthetic_process[!melanin_biosynthetic_process %in% row.names(dat)]

## testing
## 1 We create a variable with our genes' NCBI IDs
proteins <- NM_synthesis_genes_expressed
## 2 Now we map our protein IDs
proteins_mapped <- rba_string_map_ids(ids = proteins,
                                      species = 9606)

# 3 Map your IDs to STRING IDs
int_net <- rba_string_interactions_network(ids = proteins_mapped$stringId,
                                           species = 9606,
                                           required_score = 500)

# 4 plot
graph_1 <- rba_string_network_image(ids = proteins_mapped$stringId,
                                    image_format = "image",
                                    species = 9606,
                                    save_image = TRUE,
                                    required_score = 500,
                                    network_flavor = "confidence")

## Run
# Iterate through list of genes to id interactors with NM members
for (l1 in 1:length(NM_reg_list)){
  tmp <- NM_reg_list[[l1]]
  gene_compare <- tmp$Gene
  res_1 <- list()
  res_2 <- list()
  res_nm <- list()
  res_mbp <- list()
  for (i in 1:length(gene_compare)){
    # melanin_biosynthetic_process
    proteins_plus <- c(melanin_biosynthetic_process,gene_compare[i])
    proteins_mapped <- rba_string_map_ids(ids = proteins_plus,
                                          species = 9606)
    int_net <- rba_string_interactions_network(ids = proteins_mapped$stringId,
                                               species = 9606,
                                               required_score = 500)
    if(any(int_net$preferredName_A == gene_compare[i])){
      res_1[[i]] <- int_net
      res_mbp[[i]] <- sum(int_net$score[int_net$preferredName_A == gene_compare[i] | int_net$preferredName_B == gene_compare[i]])
      print(gene_compare[i])
      print(any(int_net$preferredName_A == gene_compare[i]))
    }
    
    # neuromelanin members
    proteins_plus <- c(NM_synthesis_genes_expressed,gene_compare[i])
    proteins_mapped <- rba_string_map_ids(ids = proteins_plus,
                                          species = 9606)
    int_net <- rba_string_interactions_network(ids = proteins_mapped$stringId,
                                               species = 9606,
                                               required_score = 500)
    if(any(int_net$preferredName_A == gene_compare[i])){
      res_2[[i]] <- int_net
      res_nm[[i]] <- sum(int_net$score[int_net$preferredName_A == gene_compare[i] | int_net$preferredName_B == gene_compare[i]])
      print(gene_compare[i])
      print(any(int_net$preferredName_A == gene_compare[i]))
    }
    Sys.sleep(0.1)
  }
  res_mbp[sapply(res_mbp, is.null)] <- NA
  res_nm[sapply(res_nm, is.null)] <- NA
  NM_reg_list[[l1]]$PP_MelaninBiosyntheticProcess <- unlist(lapply(res_1,is.null))
  NM_reg_list[[l1]]$PP_NeuromelaninMember <- unlist(lapply(res_2,is.null))
  
}

lapply(NM_reg_list,function(x) x["TRYP1",])


###############################################
## 2. Identify DEGs associated with Disease ###
###############################################
# 1) Annotate Voom DEGs & Find Unique DEGs accross contrasts
# voom_contrasts_interest <- c(114) # contrasts of interest eg. Diagnosis stage SNV TH
# voom_contrasts_group <- c(104:115) # contrasts of same group as interest eg. Diagnosis stage 
# 
# voom_contrasts_interest <- c(104,106,109:115) # contrasts of interest - Diagnosis stage in subregions
# voom_contrasts_group <- voom_contrasts_interest # contrasts of same group as interest - Diagnosis stage in subregions

voom_contrasts_interest <- c(28,30,33:37,39) # contrasts of interest - Diagnosis in subregions
voom_contrasts_group <- voom_contrasts_interest # contrasts of same group as interest - Diagnosis in subregions

# voom_contrasts_interest <- c(136:139,143,144) # contrasts of interest - CTR - ILBD subregions
# voom_contrasts_group <- voom_contrasts_interest # contrasts of same group as interest - CTR - ILBD subregions

# get contrasts of interest
cont_refine <- contl[contl$model.number %in% voom_contrasts_interest, ]

# annotate each contrast of interest 
cont_list <- list()
cont_uniq_list <- list()
for (z1 in 1:nrow(cont_refine)){
  # 1) annotate Voom DEGs with c-side, published cell-type DEGs, Genetic Association and published DEGs in PD
  # get matching voom and cside analysis
  c_num <- contc[unlist(lapply(strsplit(contc$limma_voom_model_compare,","), function(x) any(x %in% cont_refine$model.number[z1]))),]
  
  # get the DEG lists
  c1 <- cside_res[[c_num$model.number]]
  v1 <- voom_res[[cont_refine$model.number[z1]]]
  
  # reduce voom results to sig genes (adj p < 0.05)
  v1$Gene <- row.names(v1)
  v1 <- v1[v1$adj.P.Val < 0.05,]
  
  # reduce c-sides results to sig genes (p <0.05) and remove empty
  c1 <- lapply(c1, function(x) x[x$p_val_best < 0.05,])
  c1 <- c1[lapply(c1,nrow) > 0]
  #c1 <- lapply(c1, function(x) x[row.names(x) %in% row.names(v1),])
  
  # add cell-type to c-side results and remove zero rows
  list_of_dfs <- c1
  for (i in 1:length(list_of_dfs)){
    list_of_dfs[[i]]$cell_type <- names(list_of_dfs[i])
    list_of_dfs[[i]]$Gene <- row.names(list_of_dfs[[i]])
  }

  # combinine all data frames row-wise
  combined_df <- do.call(rbind, list_of_dfs)
  colnames(combined_df) <- paste0("cside_",colnames(combined_df))
  
  # merge with voom results
  tmp <- merge(v1,combined_df, by.x ="Gene",by.y = "cside_Gene", all.x = T)
  
  # calculate Stouffer p
  tmp$stouffer_p <- NA
  for (z in 1:nrow(tmp)){
    if(!is.na(tmp$cside_p_val_best[z])){
      tmp$stouffer_p[z] <- combine_pvalues_stouffer(c(tmp$adj.P.Val[z],tmp$cside_p_val_best[z]))
    } else {
    tmp$stouffer_p[z] <- combine_pvalues_stouffer(c(tmp$adj.P.Val[z]))
    }
  }
  combined_df <- tmp
  
  # Cross check published data cell-specific
  # create dataframe from published genesets
  deg_df <- do.call(rbind, deg_list)
  colnames(deg_df) <- paste0("Published_PD_ctDEG_evidence_",colnames(deg_df))
  # merge
  combined_df <- merge(combined_df,deg_df,by.x="Gene",by.y="Published_PD_ctDEG_evidence_Gene", all.x = T)

  # Cross check genetic associations
  gene_pd <- unlist(lapply(genetic_list, function(x) x$Gene))
  combined_df$Genetic_Association <- combined_df$Gene %in% gene_pd
  
  # Cross check publsihed data DEGs PD
  combined_df <- merge(combined_df,g4pd_deg, by.x="Gene",by.y ="Published_PD_DEG_evidence_Gene",all.x = T)
  
  # order by stouffer p and then voom p
  combined_df <- combined_df[order(combined_df$stouffer_p),]
  
  # assign to list
  cont_list[[z1]]<- combined_df
  
  # 2) Identify Voom DEGs specific to Brain region
  # retrieve contrasts of group
  cont_refine_grp <- contl[contl$model.number %in% voom_contrasts_group, ]
  
  # remove contrasts from contrast interest (both Full ROI and TH)
  cont_refine_grp <- cont_refine_grp[!cont_refine_grp$roi %in% cont_refine_grp$roi[z1],]
  
  # get contrasts with sig DEGs
  cont_refine_grp <- cont_refine_grp[cont_refine_grp$DEG_adj.p0.05 > 0,]
  
  # get the DEG lists
  r1 <- list()
  for (i in 1:nrow(cont_refine_grp)){
    r1[[i]] <- voom_res[[cont_refine_grp$model.number[i]]]
  }
  grp_deg_df <- do.call("rbind",r1)
  
  # subset for sig genes
  grp_deg_df <- grp_deg_df[grp_deg_df$adj.P.Val < 0.05,]
  
  # remove genes from contrast of interest to obtain unique DEG's
  combined_uniq_df <- combined_df[!combined_df$Gene %in% row.names(grp_deg_df),]
  
  # assign to list
  cont_uniq_list[[z1]] <- combined_uniq_df
} 

# add names to contrasts
names_1 <- paste0(cont_refine$description, "(",cont_refine$segment,")")
names(cont_list) <- names_1
names(cont_uniq_list) <- names_1

# # write to excell
# write.xlsx(cont_list, file = "Diagnosis_stage_LIMMA_Voom_annotated.xlsx", row.names = FALSE)
# write.xlsx(cont_uniq_list, file = "Diagnosis_stage_unique_LIMMA_Voom_annotated.xlsx", row.names = FALSE)

# write.xlsx(cont_list, file = "CTR.Else_LIMMA_Voom_annotated.xlsx", row.names = FALSE)
# write.xlsx(cont_uniq_list, file = "CTR.Else_unique_LIMMA_Voom_annotated.xlsx", row.names = FALSE)

# write.xlsx(cont_list, file = "CTR.ILBD_LIMMA_Voom_annotated.xlsx", row.names = FALSE)
# write.xlsx(cont_uniq_list, file = "CTR.ILBD_unique_LIMMA_Voom_annotated.xlsx", row.names = FALSE)






### Not Yet Implemented (DOWN) ####





# Identify ontologies
# count numbers of DEGs
lapply(cont_list,function(x) length(unique(x$Gene)))
lapply(cont_uniq_list,function(x) length(unique(x$Gene)))

# explore a contrast of interest
exp_cont <- as.data.frame(cont_uniq_list[[4]])
table(exp_cont$cside_cell_type)
unique(exp_cont$Gene[which(exp_cont$cside_cell_type != "<NA>")])

# get unique
tmp <- exp_cont
tmp2 <- tmp[!duplicated(tmp$Gene),]

# upreg list
upreg <- tmp2$Gene[tmp2$logFC > 0]
length(upreg)

# downreg list
downreg <- tmp2$Gene[tmp2$logFC < 0]
length(downreg)


# ontology
# select databases
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")

# upreg top 3 passing adj p-val < 0.05
enriched <- enrichr(upreg, dbs)
lapply(enriched, function(x) {
  x <- x[x$Adjusted.P.value < 0.05,]
  x[1:3,]
})

# upreg top 3 passing adj p-val < 0.05
enriched <- enrichr(downreg, dbs)
lapply(enriched, function(x) {
  x <- x[x$Adjusted.P.value < 0.05,]
  x[1:3,]
})

genes <- table(unlist(strsplit(enriched[["GO_Biological_Process_2015"]]$Genes[1:3],";")))
m_gene <- max(genes)

exp_cont[exp_cont$Gene %in% names(genes)[genes == m_gene],]



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















