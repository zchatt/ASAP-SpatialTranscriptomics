# Review of Normalization methods and their effects on GeoMx data

# Libraries
library(affy)
library(cowplot)
library(dplyr)
library(GeoMxWorkflows)
library(GeomxTools)
library(Giotto)
library(GiottoData)
library(ggforce)
library(ggnewscale)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(knitr)
library(matrixStats)
library(NanoStringNCTools)
library(plyr)
library(preprocessCore)
library(readxl)
library(reshape2)
library(Seurat)
library(scales)
library(sctransform)
library(tidyr)
library(WGCNA)
source("/Users/zacc/github_repo/Giotto/R/general_help.R")
source("/Users/zacc/github_repo/Giotto/R/utilities.R")


############################################################################################
#### Inputs
############################################################################################
run_name = "V52Y16-079-A1"
results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/visium/V52Y16-079-A1'
out_dir = "/Users/zacc/USyd/spatial_transcriptomics/data/SANPIN_VisiumFFPE_Cytassist_results/221011/V52Y16-079-A1/VISIUM/V52Y16-079-A1/outs"
analysis_dir <- results_folder

# source
source("/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/convenience/giotto_env.R")
source("/Users/zacc/github_repo/Giotto/R/general_help.R")
source("/Users/zacc/github_repo/Giotto/R/utilities.R")

############################################################################################
#### Run
############################################################################################

setwd(analysis_dir)

## load saved giotto image from "quality_control_vis.R"
visium_brain = loadGiotto(path_to_folder = file.path(results_folder, 'giotto_qc'), python_path = my_python_path)

# define raw data matrix
raw <- as.matrix(get_expression_values(visium_brain, output = "matrix",values = "raw"))

## Normalization method #1 - quantile normalization
norm.quantile <- normalize.quantiles(raw)
dimnames(norm.quantile) <- dimnames(raw)


## Normalization method #2 & #3
# SCTransform is a single-cell normalization method aimed at removing the influence of technical effects
# Also, standard log normalization for comparison
# create seurat object

seurat_obj <- Load10X_Spatial(out_dir)
seurat_obj <- SCTransform(seurat_obj , assay = "Spatial", verbose = FALSE)
seurat_obj <- NormalizeData(seurat_obj, assay = "Spatial", verbose = FALSE)

# get df
sct_df <- as.data.frame(seurat_obj@assays$SCT@counts)
log_df <- as.data.frame(seurat_obj@assays$Spatial@counts)

## Normalization method #4
# implemented as part of Giotto represents a standard procedure for normalizing for total lib size and log transformation
# Ensure the Giotto environment is installed and python path name is correct. Refer to "ASAP-SpatialTranscriptomics/convenience/giotto_env.R"
vis_tmp_df <- normalizeGiotto(gobject = visium_brain, scalefactor = 6000, verbose = T)
libsize_df <- as.matrix(get_expression_values(vis_tmp_df, output = "matrix",values = "normalized"))


### X) Comparison of normalization techniques
# get common genes and spots
common_genes <- Reduce(intersect, list(row.names(raw),
                       row.names(norm.quantile),
                       row.names(sct_df),
                       row.names(libsize_df),
                       row.names(log_df)))

common_spots <- Reduce(intersect, list(colnames(raw),
                                       colnames(norm.quantile),
                                       colnames(sct_df),
                                       colnames(libsize_df),
                                       colnames(log_df)))


# Create a list of your data frames
data_frames <- list(
  Raw = raw,
  Quantile_Norm = norm.quantile,
  SCTransform = sct_df,
  Lib_Size = libsize_df,
  Log_Norm = log_df 
)

# subset for common spots and genes
data_frames <- lapply(data_frames,function(x) x[common_genes,common_spots])
lapply(data_frames,dim)

# list all data frames
list_of_norm_frames <- data_frames
norm_names <- c("Raw","Quantile","SCTransform","Lib Size","Log_Norm")

# define dfs for plotting
ks_plot <- matrix(0,1,2)
c_plot <- matrix(0,1,2)
var_plot <- matrix(0,1,2)
pair_to_compare <- sample(1:c(ncol(data_frames[[1]]) * ncol(data_frames[[1]])), 1000, replace = FALSE) # take random set of 1000 pairs to use

for(z in 1:length(list_of_norm_frames)) {
  df_test <- as.data.frame(list_of_norm_frames[z])
  # Assess normalization methods
  tmp_grid <- expand.grid(colnames(df_test),colnames(df_test))
  tmp_grid <- tmp_grid[pair_to_compare ,] # subset for pairs to compare eg. 1000 pairs
  ks_list <- list()
  cor_list <- list()
  L <- nrow(tmp_grid)
  for(i in 1:L) {
      print(i)
      print(z)
      x <- df_test[,tmp_grid$Var1[i]] + 0.5 # add 0.5 for log trans
      y <- df_test[,tmp_grid$Var2[i]] + 0.5 # add 0.5 for log trans
       # (i) Similarity of data distributions by Kolmogorov-Smirnov Test 
      ks_list[[i]] <- ks.test(x,y)
      # (ii) Similarity of average gene expression, deviation of MA plot from y = 0
      cor_list[[i]] <- cor.test(rowMeans(log2(cbind(x,y))), c(log2(x) - log2(y)))[4]
  }
  
  ks_pval <- -log10(unlist(lapply(ks_list,function(x) as.numeric(x[2]))) + 10^-18)
  cor_val <- unlist(lapply(cor_list,function(x) as.numeric(x)))

  # calculate the row-wise variance of stable reference genes of the human brain.
  #stable_ref_genes <- c("UBE2D2", "CYC1", "RPL13") # from https://www.nature.com/articles/srep37116
  stable_ref_genes <- c("UBE2D2", "CYC1", "RPL13","TOP1","PPIA","PUM1","ACTB","UBC","B2M",
                        "ATP5B","TBP","GAPDH","EIF4A2") # using additional genes from https://www.nature.com/articles/srep37116
  
  tmp <- as.data.frame(scale(df_test[stable_ref_genes,]))
  row_var_vals <- apply(tmp,1,function(x) var(x,na.rm =T))
  row_var_vals <- row_var_vals[!is.na(row_var_vals)]
  
  
  apply(raw[row.names(raw) %in% stable_ref_genes,],1,function(x) var(x,na.rm =T))
  apply(sct_df[row.names(sct_df) %in% stable_ref_genes,],1,function(x) var(x,na.rm =T))
  
  
  # make plot dfs
  ks_plot <- rbind(ks_plot,cbind(ks_pval,rep(norm_names[z],L)))
  c_plot <- rbind(c_plot,cbind(cor_val,rep(norm_names[z],L)))
  var_plot <- rbind(var_plot,cbind(row_var_vals,rep(norm_names[z],length(stable_ref_genes))))
}

colnames(ks_plot) <- c("KS","Data")
colnames(c_plot) <- c("MAcorr","Data")
colnames(var_plot) <- c("VarRefGenes","Data")
ks_plot <- as.data.frame(ks_plot[-1,])
c_plot <- as.data.frame(c_plot[-1,])
var_plot <- as.data.frame(var_plot[-1,])

ks_plot$KS <- as.numeric(ks_plot$KS)
c_plot$MAcorr <- abs(as.numeric(c_plot$MAcorr))
c_plot <- c_plot[complete.cases(c_plot),]
var_plot$VarRefGenes <- as.numeric(var_plot$VarRefGenes)

# plots
p1 <- ggplot(ks_plot, aes(x = reorder(Data, -KS), y = KS)) +
      geom_boxplot(notch = TRUE, fill = "lightgray") +
      stat_summary(fun = mean, geom = "point",
                   shape = 18, size = 2.5, color = "#FC4E07") + 
      ylab("K-S Test -Log10(p.val)") + xlab("Data Norm Method")

p2 <- ggplot(c_plot, aes(x = reorder(Data, -MAcorr), y = MAcorr)) +
  geom_boxplot(notch = TRUE, fill = "lightgray") +
  stat_summary(fun = mean, geom = "point",
               shape = 18, size = 2.5, color = "#FC4E07") + 
  ylab("Abs(MA plot correlation)") + xlab("Data Norm Method")

p3 <- ggplot(var_plot, aes(x = reorder(Data, -VarRefGenes), y = VarRefGenes)) +
  geom_boxplot(notch = FALSE, fill = "lightgray",) +
  stat_summary(fun = mean, geom = "point",
               shape = 18, size = 2.5, color = "#FC4E07") + 
  ylab("Variance of Ref Genes") + xlab("Data Norm Method")

ggsave("normrev_3_Kolmogorov.Smirnov.test.png",p1,device = "png")
ggsave("normrev_4_MA.correlation.png",p2,device = "png")
ggsave("normrev_5_RefGenes.Variance.png",p3,device = "png")

# density distributions of all samples x genes
# Create a list to store individual density plots
density_plots <- list()
for(z in 1:length(list_of_norm_frames)) {
  print(z)
  df_test <- as.data.frame(list_of_norm_frames[z])
  df_test <- df_test[,sample(1:c(ncol(df_test)), 100, replace = FALSE)]  # sample 100 spots
  # Convert the dataframe to a long format
  long_df <- pivot_longer(df_test, everything())
  long_df$value <- long_df$value + 0.5 # adjust for log scaling
  
  # Plot the density distributions using facets
  density_plots[[z]] <- ggplot(long_df, aes(x = value, color = name)) +
    geom_density() + 
    labs(title = paste0("Method = ",norm_names[z]),
         x = "Value",
         y = "Density") +
    theme_bw() + theme(legend.position = "none") + coord_trans(x="log2")
}

p4 <- ggarrange(
  density_plots[[1]],
  density_plots[[2]],
  density_plots[[3]],
  density_plots[[4]],
  density_plots[[5]],
  ncol=2, nrow=3)

ggsave("normrev_6_density.dist.png",p4,device = "png")

