# Differential Expression Analysis
library(edgeR)
library(fgsea)
library(reactome.db)
library(metafor)
library(openxlsx)
library(viridis)
library(ggpubr)

############################################################################################
#### Inputs
############################################################################################
run_name = "geomx_oct2023"
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
rdata <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata"
contrast_path <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/contrasts_matrix.xlsx"

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
for (val in 136:nrow(cont_dat)){
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
  res[[val]] <- tmp[tmp$P.Value < 0.05,]

}


# i. name each item by description
voom_res <- res
names_1 <- paste0(cont_dat$description, "(",cont_dat$segment,")")
names(voom_res) <- names_1
for (i in 1:length(voom_res)){
  voom_res[[i]]$Contrast <- names(voom_res)[i]
}

# ## SAVE results list
# save(voom_res,cont_dat,file="/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/voom_deg.rds")


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
# colour palettes
Diagnosis_col = c("grey","#00AFBB", "#E7B800","red")
age_group_col = viridis(4)
ROI_col = viridis(7)

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
  print(stat.test)

  # add p-values to plot and save the result
  bxp <- bxp + stat_pvalue_manual(stat.test,
                                  y.position = max(data_table$y) * 1.4, step.increase = 0.1,
                                  label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)

  #return object
  return(bxp)
}

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



# Violin plot function
violin_plot_function_p <- function(data_table, x_variable, y_variable, x_lab, y_lab, colour_palette) {
  # make violin plot
  bxp <- ggviolin(
    data_table, x = x_variable, y = y_variable, 
    fill = x_variable, palette = colour_palette) + theme(legend.position = "none")
  
  # perform t-test
  data_table$x <- data_table[, x_variable]
  data_table$y <- data_table[, y_variable]
  stat.test <- data_table %>% t_test(y ~ x) %>% add_significance("p.adj")
  stat.test <- stat.test[stat.test$p.adj < 0.01,]
  print(stat.test)
  
  # add p-values to plot and save the result
  bxp <- bxp + stat_pvalue_manual(stat.test,
                                  y.position = max(data_table$y) * 1.4, step.increase = 0.1,
                                  label = "p.adj.signif") + ylab(y_lab) + xlab(x_lab)
  
  #return object
  return(bxp)
}


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
###### FGSEA
############################################################################################
# load CSIDE results; generated from ctDEG_gx, list of DEG's should be named "voom_res"
#load("voom_deg.rds")

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
msigdbr_df<- all_gene_sets %>%
  dplyr::filter(gs_cat == "H")
msigdbr_list_H = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

#msigdbr_list_reactome <- reactomePathways(row.names(exp_dat))

# # run per DEG list
# res_c2<- list()
# res_c5<- list()
# res_c8<- list()
# res_h<- list()
# #res_reactome<- list()

for (i in 96:length(voom_res)){
  print(i)
  
  # reduce to sig DEG's
  x <- voom_res[[i]]
  x <-x[x$adj.P.Val < 0.05,]
  sig_gene_list <- x[order(x$adj.P.Val),]
  
  if ( length(sig_gene_list$logFC) > 0) {
    # select res from list, perform FGSEA and plot results
    res <- as.data.frame(sig_gene_list)
    gene_stats <- res[,grep("logFC",colnames(res))]
    names(gene_stats) <- row.names(res)
    
    # FGSEA - C2
    fgseaRes <- fgsea(pathways = msigdbr_list_c2,stats= gene_stats, minSize  = 15, maxSize  = 500)
    res_c2[[i]] <- fgseaRes[padj < 0.05,]
  
    # FGSEA - C5
    fgseaRes <- fgsea(pathways = msigdbr_list_c5,stats= gene_stats, minSize  = 15, maxSize  = 500)
    res_c5[[i]] <- fgseaRes[padj < 0.05,]
    
    # FGSEA - C8
    fgseaRes <- fgsea(pathways = msigdbr_list_c8,stats= gene_stats, minSize  = 15, maxSize  = 500)
    res_c8[[i]] <- fgseaRes[padj < 0.05,]
    
    # FGSEA - H
    fgseaRes <- fgsea(pathways = msigdbr_list_H,stats= gene_stats, minSize  = 15, maxSize  = 500)
    res_h[[i]] <- fgseaRes[padj < 0.05,]
    
  }
    # # FGSEA - reactome
    # fgseaRes <- fgsea(pathways = msigdbr_list_reactome,stats= gene_stats, minSize  = 15, maxSize  = 500)
    # res_reactome[[i]] <- fgseaRes[padj < 0.05,]

}

# save fgsea results
# save(res_c2,res_c5,res_c8,res_h,file="voom_deg_fgsea.rds")

# view results
n=75

head(voom_res[[n]])

tmp <- as.data.frame(res_c2[[n]])
tmp <- tmp[tmp$padj < 0.05,]
head(tmp[order(tmp$padj),])



############################################################################################
###### Merge VOOM and FGSEA results and prioritise features
############################################################################################

### 1. summarise DEG analysis
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

# read in contrast matrix
cont_dat <- read_excel(contrast_path, sheet = "LIMMA_Voom")
cont_dat <- cont_dat[cont_dat$notes == "run",]

# add model.number to each LIMMA Voom results
list_of_dfs <- voom_res
list_of_dfs <- lapply(seq_along(list_of_dfs), function(i) {
  list_of_dfs[[i]]$model.number <- i
  list_of_dfs[[i]]$Gene <- row.names(list_of_dfs[[i]])
  return(list_of_dfs[[i]])
})

# Combining all data frames row-wise
combined_df <- do.call(rbind, list_of_dfs)

# merge modelling information
combined_df <- merge(combined_df,cont_dat,by="model.number")

# remove models not wanted
combined_df <- combined_df[!combined_df$model.number %in% c(1:19),] # Brainregion models with all samples, not just controls

# sort by p-value
combined_df <- combined_df[combined_df$adj.P.Val < 0.05,]
combined_df <- combined_df[order(combined_df$adj.P.Val),]

# define each genes most sig contrast for each concept
gene_id <- unique(combined_df$Gene)
res_agg <- list()
for (i in 1:length(gene_id)){
  print(i)
  # get top p and all analysis descriptions
  tmp <- combined_df[combined_df$Gene == gene_id[i] & combined_df$concept_bin == "Brainregion",]
  Brainregion_contrast <- paste(paste0(tmp$description," [",tmp$segment,"]"), collapse = " | ")
  Brainregion_contrast_top_p <- tmp$adj.P.Val[1]
  
  tmp <- combined_df[combined_df$Gene == gene_id[i] & combined_df$concept_bin == "Disease",]
  Disease_contrast <- paste(paste0(tmp$description," [",tmp$segment,"]"), collapse = " | ")
  Disease_contrast_top_p <- tmp$adj.P.Val[1]
  
  tmp <- combined_df[combined_df$Gene == gene_id[i] & combined_df$concept_bin == "Neuromelanin",]
  Neuromelanin_contrast <- paste(paste0(tmp$description," [",tmp$segment,"]"), collapse = " | ")
  Neuromelanin_contrast_top_p <- tmp$adj.P.Val[1]
  
  # create stouffer p-value
  p_values <- c(Brainregion_contrast_top_p,Disease_contrast_top_p, Neuromelanin_contrast_top_p)
  p_values <- p_values[!is.na(p_values)]
  concepts_stouffer_p <- combine_pvalues_stouffer(p_values)
  
  res_agg[[i]] <- cbind(Gene = gene_id[i],Brainregion_contrast, Brainregion_contrast_top_p,
                        Disease_contrast, Disease_contrast_top_p,
                        Neuromelanin_contrast, Neuromelanin_contrast_top_p, concepts_stouffer_p)
}

# Combining all data frames row-wise
res_agg_df <- as.data.frame(do.call(rbind, res_agg))














### 2. summarise FGSEA 

# add model.number to each LIMMA Voom results
list_of_dfs <- res_c2
list_of_dfs <- lapply(seq_along(list_of_dfs), function(i) {
  list_of_dfs[[i]]$model.number <- i
  return(list_of_dfs[[i]])
})

gene_id[[i]]
res_c2

lapply(res_c2, function(i) {
  
  x[x$ES < 0,][1:5]
  x[x$ES > 0,][1:5]
  list_of_dfs[[i]]$model.number <- i
  list_of_dfs[[i]]$Gene <- row.names(list_of_dfs[[i]])
  return(list_of_dfs[[i]])
})