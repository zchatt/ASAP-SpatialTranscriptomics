# human evaluation of CHCHD2/10 and SNCA expression within SNV & SND DA neurons of controls and PD by GeoMx Spatial Transcriptomics

library(ggplot2)
library(dplyr)
library(viridis)
library(data.table)
library(ggpubr)
library(rstatix)
library(spacexr)
library(Seurat)
library(edgeR)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(dplyr)


#----------- Input ----------- 

analysis_dir = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis"
setwd(analysis_dir)

load("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/geomx_oct2023_seurat.Rdata") # ST - GeoMx thresheld data


#----------- Cohort stats ----------- 
tmp <- meta_dat[!duplicated(meta_dat$Brainbank_ID),]

#Sex
tmp %>%
  group_by(Diagnosis,Sex) %>%
  dplyr::summarise(count = n_distinct(Brainbank_ID))


meta_dat_old <- read_excel("/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_oct2023/analysis/targeted_genelist_geomx_011223.xlsx")
tmp <- meta_dat_old[!duplicated(meta_dat_old$Brainbank_ID),]

#Sex
tmp %>%
  group_by(Diagnosis,Sex) %>%
  dplyr::summarise(count = n_distinct(Brainbank_ID))


#----------- DEG ----------- 
# extract data matrices - counts for limma
exp_dat <- as.matrix(gxdat_s@assays$RNA$counts)
meta_dat <- as.data.frame(gxdat_s@meta.data)
colnames(meta_dat) <- make.names(colnames(meta_dat))

genes_interest <- c("SNCA","CHCHD2","CHCHD10")

# define df and meta data
meta_dat <- meta_dat[meta_dat$segment == "TH" & # DA neurons
                       meta_dat$ROI %in% c("SNV","SND") & # SNV and SND restions
                       meta_dat$Diagnosis %in% c("CTR","ILBD","ePD","lPD"),] # Dx
exp_dat <- exp_dat[,row.names(meta_dat)]

# format covariates
meta_dat$Sex <- as.factor(meta_dat$Sex)
meta_dat$Age  <- as.numeric(meta_dat$Age)
meta_dat$DV200  <- as.numeric(meta_dat$DV200)
meta_dat$PMD.hs <- as.numeric(meta_dat$PMD.hs)
meta_dat$Brainbank_ID <- as.factor(meta_dat$Brainbank_ID)
meta_dat$scan.name <- as.factor(meta_dat$scan.name)

## 1)  limma voom - SNV v SND within Diagnosis
contrast_var <- unique(meta_dat$Diagnosis)
meta_dat$condt <- meta_dat$ROI == "SNV"
  
res <- list()
for (val in 1:length(contrast_var)){  
  print(contrast_var[val])
  
  # subset for var of interest
  meta_dat2 <- meta_dat[meta_dat$Diagnosis == contrast_var[val],]
  exp_dat2 <- exp_dat[,row.names(meta_dat2)]
  meta_dat2 <- droplevels(meta_dat2)
  
  dge <- DGEList(exp_dat2, group = meta_dat2$condt)
  dge <- calcNormFactors(dge)
  design <- model.matrix( ~ condt + Brainbank_ID, data=meta_dat2)
  vm <- voom(dge, design = design, plot = TRUE)
  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)
  tt <- topTable(fit, n = Inf, adjust.method = "BH")
  
  tt$contrast <- contrast_var[val]
  res[[val]] <- tt[genes_interest,]
}

# format results
res <- lapply(res, function(x) {
  x <- x[,-grep("Brainbank_ID",colnames(x))]
  x$Gene <- row.names(x)
  row.names(x) <- NULL
  return(x)
})
stat.summary <- do.call(rbind,res)
df <- as.data.frame(apply(stat.summary,2,as.character))

df$stars <- cut(as.numeric(df$adj.P.Val),
                breaks = c(-Inf, 0.001, 0.01, 0.05,0.1, Inf),
                labels = c("***", "**", "*", ".",""),
                right = FALSE)

# write summary to file
write.table(df, file="deg_snvVsnd.dx_THpos.stat.summary.txt", sep="\t",row.names = F, quote = F)



## 2) limma voom - SNV between ILBD/ ePD/ lPD and Control
contrast_var <- unique(meta_dat$Diagnosis)
contrast_var <- contrast_var[!contrast_var == "CTR"]
genes_interest <- c("SNCA","CHCHD2","CHCHD10")

res <- list()
for (val in 1:length(contrast_var)){  
  print(contrast_var[val])
  
  # subset for var of interest
  meta_dat2 <- meta_dat[meta_dat$Diagnosis %in% c("CTR",contrast_var[val]) & meta_dat$ROI == "SND",]
  exp_dat2 <- exp_dat[,row.names(meta_dat2)]
  meta_dat2$condt <- meta_dat2$Diagnosis == "CTR"
  meta_dat2 <- droplevels(meta_dat2)
  
  dge <- DGEList(exp_dat2, group = meta_dat2$condt)
  dge <- calcNormFactors(dge)
  design <- model.matrix( ~ condt + Age + PMD.hs + Sex, data=meta_dat2)
  vm <- voom(dge, design = design, plot = TRUE)
  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)
  tt <- topTable(fit, n = Inf, adjust.method = "BH")
  
  tt$contrast <- contrast_var[val]
  res[[val]] <- tt[genes_interest,]
}

# format results
res <- lapply(res, function(x) {
  x$Gene <- row.names(x)
  row.names(x) <- NULL
  return(x)
})
stat.summary <- do.call(rbind,res)
df <- as.data.frame(apply(stat.summary,2,as.character))

df$stars <- cut(as.numeric(df$adj.P.Val),
                breaks = c(-Inf, 0.001, 0.01, 0.05,0.1, Inf),
                labels = c("***", "**", "*", ".",""),
                right = FALSE)

# write summary to file
#write.table(df, file="deg_snv.ctrVdx_THpos.stat.summary.txt", sep="\t",row.names = F, quote = F)
write.table(df, file="deg_snd.ctrVdx_THpos.stat.summary.txt", sep="\t",row.names = F, quote = F)



#------------ mean and sd ----------
# extract data matrices - normalized counts for plotting
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)

# define genes to plot
genes_to_plot <- c("SNCA","CHCHD2","CHCHD10")

# define metadata
meta_dat1 <- meta_dat[meta_dat$ROI %in% c("SND") &
                        meta_dat$segment == "TH" &
                        meta_dat$Diagnosis %in% c("CTR","ILBD","ePD","lPD"),]
exp_dat <- exp_dat[genes_to_plot,row.names(meta_dat1)]

# combine data
data_table <- cbind(meta_dat1,t(exp_dat))

# mean and sd
data_table %>% group_by(Diagnosis) %>%
  summarise(mean=mean(CHCHD2), sd=sd(CHCHD2))

data_table %>% group_by(Diagnosis) %>%
  summarise(mean=mean(CHCHD10), sd=sd(CHCHD10))

data_table %>% group_by(Diagnosis) %>%
  summarise(mean=mean(SNCA), sd=sd(SNCA))


#------------tabulation of cohort statistics ------------
# create a glm of cohort stats
tmp <- meta_dat

# Age
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(median_age = median(Age,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(iqr_age = IQR(Age,na.rm=TRUE))


#Sex
tmp %>%
  group_by(Diagnosis,Sex) %>%
  dplyr::summarise(count = n_distinct(Brainbank_ID))

tmp2 <- tmp[!duplicated(tmp$Brainbank_ID),]
chisq.test(tmp2$Sex,tmp2$Diagnosis)


# postmortem interval
tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(result = median(PMD.hs,na.rm=TRUE))

tmp %>%
  distinct(Brainbank_ID, .keep_all=TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(result = IQR(PMD.hs,na.rm=TRUE))


## compare tabulates data 
tmp <- meta_dat[!duplicated(meta_dat$Brainbank_ID),]

# Loop over diagnoses and compare to "CTR"
diagnoses <- c("ILBD", "ePD", "lPD")

# Initialize lists to store results
chi_results <- list()
t_age_results <- list()
t_pmd_results <- list()

# Function for individual comparisons
for (diag in diagnoses) {
  # Subset the data for "CTR" and the current diagnosis
  tmp_sub <- subset(tmp, Diagnosis %in% c("CTR", diag))
  
  # Chi-square test for "Sex"
  sex_table <- table(tmp_sub$Sex, tmp_sub$Diagnosis)
  chi_results[[diag]] <- chisq.test(sex_table)
  print(paste("Chi-square test for 'Sex':", diag, "vs CTR"))
  print(chi_results[[diag]])
  
  # Test for 'Age'
  if (shapiro.test(tmp_sub$Age)$p.value > 0.05) {
    t_age_results[[diag]] <- t.test(Age ~ Diagnosis, data = tmp_sub)
  } else {
    t_age_results[[diag]] <- wilcox.test(Age ~ Diagnosis, data = tmp_sub)
  }
  print(paste("Test for 'Age':", diag, "vs CTR"))
  print(t_age_results[[diag]])
  
  # Test for 'PMD.hs'
  if (shapiro.test(tmp_sub$PMD.hs)$p.value > 0.05) {
    t_pmd_results[[diag]] <- t.test(PMD.hs ~ Diagnosis, data = tmp_sub)
  } else {
    t_pmd_results[[diag]] <- wilcox.test(PMD.hs ~ Diagnosis, data = tmp_sub)
  }
  print(paste("Test for 'PMD.hs':", diag, "vs CTR"))
  print(t_pmd_results[[diag]])
}

# Apply Bonferroni correction for Age and PMD.hs
# Since there are 3 comparisons, the correction factor is 3
bonferroni_age_p <- sapply(t_age_results, function(x) x$p.value * 3)
bonferroni_pmd_p <- sapply(t_pmd_results, function(x) x$p.value * 3)

print("Bonferroni corrected p-values for Age:")
print(bonferroni_age_p)

print("Bonferroni corrected p-values for PMD.hs:")
print(bonferroni_pmd_p)



#------------ Barplots ----------
# extract data matrices - normalized counts for plotting
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)

# define genes to plot
genes_to_plot <- c("SNCA","CHCHD2","CHCHD10")

# define metadata
meta_dat1 <- meta_dat[meta_dat$ROI %in% c("SNV","SND") &
                       meta_dat$segment == "TH" &
                       meta_dat$Diagnosis %in% c("CTR","ILBD","ePD","lPD"),]
exp_dat <- exp_dat[genes_to_plot,row.names(meta_dat1)]

# combine data
data_table <- cbind(meta_dat1,t(exp_dat))

# Define a function to calculate SEM
calculate_sem <- function(x) {
  return(sd(x) / sqrt(length(x)))
}

# Custom color mapping for Diagnosis
color_mapping <- c("SND" = "red", "SNV" = "blue")

#detach(package:plyr)
# Create a plotting function for each gene
plot_gene_expression <- function(gene) {
  # Specify the desired order of Diagnosis levels
  data_table$Diagnosis[data_table$Diagnosis == "ILBD"] <- "ILB"
  data_table$Diagnosis <- factor(data_table$Diagnosis, levels = c("CTR","ILB","ePD","lPD" ))  # Adjust levels as needed
  
  # Summarize the data
  data_summary <- data_table %>%
    group_by(Diagnosis, ROI) %>%
    summarise(
      Mean_Expression = mean(!!sym(gene), na.rm = TRUE),
      SEM_Expression = calculate_sem(!!sym(gene))
    )
  
  # Plot the data
  ggplot(data_summary, aes(x = Diagnosis, y = Mean_Expression, fill = ROI)) +
    geom_bar(stat = "identity", position = position_dodge(0.9), color = "black", width = 0.5) +  # Black outline
    geom_errorbar(aes(ymin = Mean_Expression, ymax = Mean_Expression + SEM_Expression),
                  position = position_dodge(0.9), width = 0.25) +  # Only positive SEM bars
    scale_fill_manual(values = color_mapping) +  # Custom colors
    labs(title = paste(gene, "Expression"), y = "Normalized Expression", x = "") +
    theme_classic() +
    #ylim(0, max(data_summary$Mean_Expression) * 1.4)
    ylim(0, 75)
  
}

# Plot for SNCA
plot_snca <- plot_gene_expression("SNCA")

# Plot for CHCHD2
plot_chchd2 <- plot_gene_expression("CHCHD2")

# Plot for CHCHD10
plot_chchd10 <- plot_gene_expression("CHCHD10")

# arrange
arrange <- ggarrange(plotlist=list(plot_chchd2, plot_chchd10,plot_snca), nrow=2, ncol=3)
ggsave("barplots_chchd2.chchd10.snca_dx.snv.snd_normexp2.pdf", arrange,width = 8, height = 6)





### Sex
# extract data matrices - normalized counts for plotting
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)

# define genes to plot
genes_to_plot <- c("SNCA","CHCHD2","CHCHD10")

# define metadata
meta_dat1 <- meta_dat[meta_dat$ROI %in% c("SNV","SND") &
                       meta_dat$segment == "TH" &
                       meta_dat$Diagnosis %in% c("CTR"),]
exp_dat <- exp_dat[genes_to_plot,row.names(meta_dat1)]

# combine data
data_table <- cbind(meta_dat1,t(exp_dat))

# Define a function to calculate SEM
calculate_sem <- function(x) {
  return(sd(x) / sqrt(length(x)))
}

# Custom color mapping for Diagnosis
color_mapping <- c("SND" = "red", "SNV" = "blue")

#detach(package:plyr)
# Create a plotting function for each gene
plot_gene_expression <- function(gene) {
  # Specify the desired order of Diagnosis levels
  #data_table$Diagnosis <- factor(data_table$Diagnosis, levels = c("CTR","ILBD","ePD","lPD" ))  # Adjust levels as needed
  
  # Summarize the data
  data_summary <- data_table %>%
    group_by(Sex, ROI) %>%
    summarise(
      Mean_Expression = mean(!!sym(gene), na.rm = TRUE),
      SEM_Expression = calculate_sem(!!sym(gene))
    )
  
  # Plot the data
  ggplot(data_summary, aes(x = Sex, y = Mean_Expression, fill = ROI)) +
    geom_bar(stat = "identity", position = position_dodge(0.9), color = "black", width = 0.5) +  # Black outline
    geom_errorbar(aes(ymin = Mean_Expression, ymax = Mean_Expression + SEM_Expression),
                  position = position_dodge(0.9), width = 0.25) +  # Only positive SEM bars
    scale_fill_manual(values = color_mapping) +  # Custom colors
    labs(title = paste(gene, "Expression"), y = "Normalized Expression", x = "") +
    theme_classic() +
    ylim(0, max(data_summary$Mean_Expression) * 1.5)
}

# Plot for SNCA
plot_snca <- plot_gene_expression("SNCA")

# Plot for CHCHD2
plot_chchd2 <- plot_gene_expression("CHCHD2")

# Plot for CHCHD10
plot_chchd10 <- plot_gene_expression("CHCHD10")

# arrange
arrange <- ggarrange(plotlist=list(plot_chchd2, plot_chchd10,plot_snca), nrow=2, ncol=3)
#ggsave("barplots_chchd2.chchd10.snca_dx.snv.snd_normexp.pdf", arrange,width = 8, height = 6)




#------------ Correlation plots ----------
# extract data matrices - normalized counts for plotting
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)

# define genes to plot
genes_to_plot <- c("SNCA","CHCHD2","CHCHD10")

# define metadata
meta_dat1 <- meta_dat[meta_dat$ROI %in% c("SNV") &
                        meta_dat$segment == "TH" &
                        meta_dat$Diagnosis %in% c("CTR"),]
data_table <- as.data.frame(cbind(meta_dat1,t(exp_dat[genes_to_plot,row.names(meta_dat1)])))

# Create scatterplot for CHCHD2 vs SNCA and CHCHD10 vs SNCA with regression lines
s1 <- ggplot(data_table, aes(x = SNCA)) +
  
  # Scatter plot and regression line for CHCHD2 vs SNCA
  geom_point(aes(y = CHCHD2), color = "blue", alpha = 0.6, shape = 19) +
  geom_smooth(aes(y = CHCHD2), method = "lm", color = "blue", se = FALSE) +
  
  # Add R and p-value for CHCHD2 vs SNCA
  stat_cor(
    aes(y = CHCHD2,label = paste0("CHCHD2:", ..rr.label..,sep = "~`,`~", ..p.label..)), 
    label.x = 3, label.y = max(data_table$CHCHD10) + 20, color = "blue"
  ) +

  # Scatter plot and regression line for CHCHD10 vs SNCA
  geom_point(aes(y = CHCHD10), color = "darkorchid", alpha = 0.6 , shape = 19) +
  geom_smooth(aes(y = CHCHD10), method = "lm", color = "darkorchid", se = FALSE) +
  
  # Add R and p-value for CHCHD10 vs SNCA
  stat_cor(
    aes(y = CHCHD10,label = paste0("CHCHD10:", ..rr.label..,sep = "~`,`~", ..p.label..)), 
    label.x = 3, label.y = max(data_table$CHCHD10) + 9 , color = "darkorchid"
  ) +
  
  # Titles and labels
  labs(
    title = "SNV",
    x = "SNCA Expression",
    y = "CHCHD2 or CHCHD10 Norm. Expression"
  ) +
  
  # Theme adjustments for better visualization
  theme_classic()


# define metadata
meta_dat1 <- meta_dat[meta_dat$ROI %in% c("SND") &
                        meta_dat$segment == "TH" &
                        meta_dat$Diagnosis %in% c("CTR"),]
data_table <- as.data.frame(cbind(meta_dat1,t(exp_dat[genes_to_plot,row.names(meta_dat1)])))

# Create scatterplot for CHCHD2 vs SNCA and CHCHD10 vs SNCA with regression lines
s2 <- ggplot(data_table, aes(x = SNCA)) +
  
  # Scatter plot and regression line for CHCHD2 vs SNCA
  geom_point(aes(y = CHCHD2), color = "blue", alpha = 0.6, shape = 19) +
  geom_smooth(aes(y = CHCHD2), method = "lm", color = "blue", se = FALSE) +
  
  # Add R and p-value for CHCHD2 vs SNCA
  stat_cor(
    aes(y = CHCHD2,label = paste0("CHCHD2:", ..rr.label..,sep = "~`,`~", ..p.label..)), 
    label.x = 3, label.y = max(data_table$CHCHD10) + 20, color = "blue"
  ) +
  
  # Scatter plot and regression line for CHCHD10 vs SNCA
  geom_point(aes(y = CHCHD10), color = "darkorchid", alpha = 0.6 , shape = 19) +
  geom_smooth(aes(y = CHCHD10), method = "lm", color = "darkorchid", se = FALSE) +
  
  # Add R and p-value for CHCHD10 vs SNCA
  stat_cor(
    aes(y = CHCHD10,label = paste0("CHCHD10:", ..rr.label..,sep = "~`,`~", ..p.label..)), 
    label.x = 3, label.y = max(data_table$CHCHD10) + 9 , color = "darkorchid"
  ) +
  
  # Titles and labels
  labs(
    title = "SND",
    x = "SNCA Expression",
    y = "CHCHD2 or CHCHD10 Norm. Expression"
  ) +
  
  # Theme adjustments for better visualization
  theme_classic()

# arrange
arrange <- ggarrange(plotlist=list(s1,s2), nrow=2, ncol=2)
ggsave("scatter_chchd2.chchd10.snca_snv.snd.controls_normexp.pdf", arrange,width = 8, height = 6)
