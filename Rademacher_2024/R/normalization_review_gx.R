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


# # # Input
# # Locally produced data from Substantia Nigra
run_name = "CHA12467"
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/CHA12467/analysis"

# run_name <- "hu_brain_001"
# analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/hu_brain_workflow_and_counts_files/analysis"
# my_python_path = '/Users/zacc/opt/anaconda3/envs/py311/bin/python'

#########
## RUN ##
#########
setwd(analysis_dir)
# Load qc data generated from geomx_lowlevel.R, should be in "/analysis" of the name "run_name_qc.gx.Rdata"
load(paste0(run_name, "_qc.gx.Rdata"))

## Normalization method #1
# 1) Using i) quartile 3 (Q3) or ii) background normalization.
# From https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

#  manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "DetectionThreshold"
Stat_data <- data.frame(row.names = colnames(exprs(target_demoData)),
                        Segment = colnames(exprs(target_demoData)),
                        Annotation = pData(target_demoData)[, ann_of_interest],
                        Q3 = unlist(apply(exprs(target_demoData), 2, quantile, 0.75, na.rm = TRUE)),
                        NegProbe = exprs(target_demoData)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43, 0.57))

png("normrev_1_Q3vsnegGeoMean.png")
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
dev.off()

# define raw data matrix
raw <- as.matrix(exprs(target_demoData))

# Q3 norm (75th percentile) for WTA/CTA with or without custom spike-ins
target_demoData <- normalize(target_demoData ,
                             norm_method = "quant",
                             desiredQuantile = .75,
                             toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_demoData <- normalize(target_demoData ,
                             norm_method = "neg",
                             fromElt = "exprs",
                             toElt = "neg_norm")

# define Q3 normalized and background normalization
q3_norm <- assayDataElement(target_demoData, elt = "q_norm")
back_norm <- assayDataElement(target_demoData, elt = "neg_norm")

## Normalization method #2 - quantile normalization
# Recently Van Hijfte et al. (iScience. 2023) reviewed the Nanostring GeoMx Q3 normalization technique - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9800292/
# Observed that Q3 normalization does not correct for large differences in the
# signal (gene expression) to noise (neg probes) observed between samples.
# Authors assessed CPM normalization, DESeq2 normalization, gamma fit correction, and quantile normalization (id as best).
norm.quantile <- normalize.quantiles(raw)
dimnames(norm.quantile) <- dimnames(exprs(target_demoData))


## Normalization method #3
# SCTransform is a single-cell normalization method aimed at removing the influence of technical effects
# create seurat object
geomx_raw <- CreateSeuratObject(raw,
  project = "RawGeoMx",
  assay = "RNA",
  names.field = 1)

# run sctransform
geomx_raw <- SCTransform(geomx_raw)

# get df
geomx_sct_df <- as.data.frame(geomx_raw@assays$SCT@counts)


## Normalization method #4
# implemented as part of Giotto represents a standard procedure for normalizing for total lib size and log transformation
# Ensure the Giotto environment is installed and python path name is correct. Refer to "ASAP-SpatialTranscriptomics/convenience/giotto_env.R"
instrs = createGiottoInstructions(python_path = my_python_path)
# Create a Giotto object
giotto_df = createGiottoObject(expression = raw, instructions = instrs)

## normalize
giotto_df <- normalizeGiotto(gobject = giotto_df, scalefactor = 6000, verbose = T)
libsize_df <- as.matrix(get_expression_values(giotto_df, output = "matrix",values = "normalized"))

### X) Plotting counts  of normalization techniques
b1 <- boxplot(raw, col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:ncol(raw), xlab = "Segment",
        ylab = "Counts, Raw")

b2 <- boxplot(q3_norm, col = "#9EDAE5", main = "Q3 Norm Counts",
              log = "y", names = 1:ncol(raw), xlab = "Segment",
              ylab = "Counts, Q3 Normalized")

b3 <- boxplot(back_norm, col = "#9EDAE5", main = "Neg Norm Counts",
              log = "y", names = 1:ncol(raw), xlab = "Segment",
              ylab = "Counts, Neg. Normalized")

b4 <- boxplot(norm.quantile, col = "#9EDAE5", main = "Quantile Norm Counts",
              log = "y", names = 1:ncol(raw), xlab = "Segment",
              ylab = "Counts, Quantile Normalized")

b5 <- boxplot(geomx_sct_df, col = "#9EDAE5", main = "SCTransform Norm Counts",
              log = "y", names = 1:ncol(raw), xlab = "Segment",
              ylab = "Counts, SCTransform Normalized")

b6 <- boxplot(libsize_df, col = "#9EDAE5", main = "Lib Size Norm Counts",
              log = "y", names = 1:ncol(raw), xlab = "Segment",
              ylab = "Counts, Lib Size Normalized")

# Create a list of your data frames
data_frames <- list(
  Raw = raw,
  Q3_Norm = q3_norm,
  Neg_Norm = back_norm,
  Quantile_Norm = norm.quantile,
  SCTransform = geomx_sct_df,
  Lib_Size = libsize_df
)

# Create a list to store the plots
plot_list <- list()

# Create the boxplot
for (i in seq_along(data_frames)) {
  all_data <- as.data.frame(data_frames[[i]])
  tidy_data <- all_data %>%
    pivot_longer(names(.), names_to = "Segment", values_to = "Counts")
  
  p <- ggplot(tidy_data, aes(x = Segment, y = Counts)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(title = names(data_frames)[i], x = "Segment", y = "Counts") +
    theme_minimal()
  
plot_list[[i]] <- p
}

# Arrange the plots using ggarrange
combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = 2, nrow = 3)

png("normrev_2_countsnorm.png")
combined_plot 
dev.off()

### X) Comparison of normalization techniques
# list all data frames
list_of_norm_frames <- data_frames
norm_names <- c("Raw","Q3","Backg Norm","Quantile","SCTransform","Lib Size")

# define dfs for plotting
ks_plot <- matrix(0,1,2)
c_plot <- matrix(0,1,2)
var_plot <- matrix(0,1,2)

for(z in 1:length(list_of_norm_frames)) {
  df_test <- as.data.frame(list_of_norm_frames[z])
  # Assess normalization methods
  tmp_grid <- expand.grid(colnames(df_test),colnames(df_test))
  ks_list <- list()
  cor_list <- list()
  L <- nrow(tmp_grid)
  for(i in 1:L) {
      print(i)
      print(z)
      x <- df_test[,tmp_grid$Var1[i]]
      y <- df_test[,tmp_grid$Var2[i]]
       # (i) Similarity of data distributions by Kolmogorov-Smirnov Test 
      ks_list[[i]] <- ks.test(x,y)
      # (ii) Similarity of average gene expression, deviation of MA plot from y = 0
      cor_list[[i]] <- cor.test(rowMeans(log2(cbind(x,y))), c(log2(x) - log2(y)))[4]
  }
  
  ks_pval <- -log10(unlist(lapply(ks_list,function(x) as.numeric(x[2]))) + 10^-18)
  cor_val <- unlist(lapply(cor_list,function(x) as.numeric(x)))

  # calculate the row-wise variance of stable reference genes of the human brain.
  stable_ref_genes <- c("UBE2D2", "CYC1", "RPL13")
  row_var_vals <- rowVars(as.matrix(scale(df_test[stable_ref_genes,])))
  
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
  df_test <- as.data.frame(list_of_norm_frames[z])
  # Convert the dataframe to a long format
  long_df <- pivot_longer(df_test, everything())
  long_df$value <- long_df$value + 0.5 # adust for log scaling
  
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
  density_plots[[6]],
  ncol=2, nrow=3)

ggsave("normrev_6_density.dist.png",p4,device = "png")

