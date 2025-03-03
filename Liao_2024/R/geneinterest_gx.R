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

## boxplots 
dplot <- as.matrix(gxdat_s@assays$RNA$data)[genes_interest ,row.names(meta_dat)]
dplot <- cbind(t(dplot),meta_dat[,c("Diagnosis","ROI")])
dplot <- dplot[dplot$Diagnosis == "CTR",]

# Ensure ROI is a factor for proper ordering
dplot$ROI <- factor(dplot$ROI, levels = c("SND", "SNV"))

# Convert to long format for faceting
dplot_long <- dplot %>%
  pivot_longer(cols = c(SNCA, CHCHD2, CHCHD10), 
               names_to = "Gene", 
               values_to = "Expression")

# Prepare p-values for annotation
df$adj.P.Val <- as.numeric(df$adj.P.Val)

# Prepare p-values for annotation
df_pvalues <- df %>%
  filter(contrast == "CTR") %>%
  select(Gene, adj.P.Val, stars) %>%
  mutate(label = stars)

# Merge p-values into the plotting data
dplot_long <- dplot_long %>%
  left_join(df_pvalues, by = "Gene")

# Generate the boxplot with fixed y-axis scales and custom colors
g1 <- ggplot(dplot_long, aes(x = ROI, y = Expression, fill = ROI)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~ Gene) +  # Use scales = "free_y" if each gene has a different range
  scale_fill_manual(values = c("SND" = "grey", "SNV" = "black")) +
  labs(x = "", y = "Normalized Expression Level", 
       title = "") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  geom_text(data = df_pvalues, aes(x = 1.5, y = Inf, label = label), 
            inherit.aes = FALSE, vjust = 1.5, size = 5)

g2 <- ggplot(dplot_long, aes(x = ROI, y = Expression, fill = ROI)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~ Gene, scales = "free_y") +  # Use scales = "free_y" if each gene has a different range
  scale_fill_manual(values = c("SND" = "grey", "SNV" = "black")) +
  labs(x = "", y = "Normalized Expression Level", 
       title = "") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  geom_text(data = df_pvalues, aes(x = 1.5, y = Inf, label = label), 
            inherit.aes = FALSE, vjust = 1.5, size = 5)


# arrange
arrange <- ggarrange(plotlist=list(g1,g2), nrow=2, ncol=2)
ggsave("boxplots_ctr_genes.interest_snv.snd.pdf", arrange,width = 8, height = 6)






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

tmp <- meta_dat

#Sex
tmp %>%
  group_by(Diagnosis,Sex) %>%
  dplyr::summarise(count = n_distinct(Brainbank_ID))


# Age
tmp %>%
  distinct(Brainbank_ID, .keep_all = TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(
    median_age = median(Age, na.rm = TRUE),
    Q1 = quantile(Age, 0.25, na.rm = TRUE),
    Q3 = quantile(Age, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1  # Optional: include IQR as a separate column
  )

# postmortem interval
tmp %>%
  distinct(Brainbank_ID, .keep_all = TRUE) %>%
  group_by(Diagnosis) %>%
  dplyr::summarize(
    median_age = median(PMD.hs, na.rm = TRUE),
    Q1 = quantile(PMD.hs, 0.25, na.rm = TRUE),
    Q3 = quantile(PMD.hs, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1  # Optional: include IQR as a separate column
  )

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




#------------ Review Q - OXPHOS subunit expression correlation to CHCHD2 and CHCHD10 ----------
library(msigdbr)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ComplexHeatmap)
library(pheatmap)

# ---- Load expression and metadata matrices
exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat <- as.data.frame(gxdat_s@meta.data)
colnames(meta_dat) <- make.names(colnames(meta_dat))

# Filter metadata for relevant samples
meta_dat1 <- meta_dat[meta_dat$ROI %in% c("SNV") &
                        meta_dat$segment == "TH" &
                        meta_dat$Diagnosis %in% c("CTR"), ]
exp_dat <- exp_dat[, row.names(meta_dat1)]


## ---- Get Oxphos genes
# Get the MSigDB hallmark gene sets
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")

# Filter for the oxidative phosphorylation hallmark gene set
oxphos_genes <- hallmark_gene_sets %>%
  filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  select(human_gene_symbol)

## cor-cor
# Define genes of interest
genes_of_interest <- c("SNCA", "CHCHD2", "CHCHD10")

# Retrieve the oxidative phosphorylation (OXPHOS) genes from the extracted gene set
oxphos_genes <- oxphos_genes$human_gene_symbol

# Combine genes for analysis
all_genes <- unique(c(oxphos_genes, genes_of_interest))

# Filter the expression data for selected genes
filtered_exp <- exp_dat[rownames(exp_dat) %in% all_genes, ]

# Check for missing genes
missing_genes <- setdiff(all_genes, rownames(filtered_exp))
if (length(missing_genes) > 0) {
  message("Warning: The following genes are not in the dataset: ", paste(missing_genes, collapse = ", "))
}

# Calculate correlation matrix
cor_matrix <- cor(t(filtered_exp), method = "pearson", use = "pairwise.complete.obs")


# Plot heatmap
library(pheatmap)
library(viridis)  # For magma color palette

## Plot heatmap
# Define genes of interest with specific colors
genes_of_interest <- c("SNCA", "CHCHD2", "CHCHD10")

# Create annotation data frame for highlighting
annotation_row <- data.frame(Gene = rownames(cor_matrix))
annotation_row$Category <- "Oxphos Genes"
annotation_row$Category[annotation_row$Gene == "SNCA"] <- "SNCA"
annotation_row$Category[annotation_row$Gene == "CHCHD10"] <- "CHCHD10"
annotation_row$Category[annotation_row$Gene == "CHCHD2"] <- "CHCHD2"
rownames(annotation_row) <- rownames(cor_matrix)  # Correct rownames
annotation_row <- annotation_row[, "Category", drop = FALSE]

# Define annotation colors for gene highlights
annotation_colors <- list(
  Category = c(
    SNCA = "red",
    CHCHD2 = "green1",
    CHCHD10 = "green4",
    "Oxphos Genes" = "grey"
  )
)

# Define correlation color palette
color_palette <- viridis::magma(50)

# Perform hierarchical clustering for rows and columns
dist_rows <- dist(cor_matrix)  # Distance matrix for rows (genes)
hclust_rows <- hclust(dist_rows, method = "ward.D2")  # Clustering for rows

dist_cols <- dist(t(cor_matrix))  # Distance matrix for columns (samples)
hclust_cols <- hclust(dist_cols, method = "ward.D2")  # Clustering for columns

# Cut the dendrogram to obtain k=2 clusters for rows and columns
cluster_rows <- cutree(hclust_rows, k = 2)
cluster_cols <- cutree(hclust_cols, k = 2)

# Add row cluster information to annotation
annotation_row$Cluster <- as.factor(cluster_rows)

# Define cluster colors for rows and columns
annotation_colors$Cluster <- c("1" = "purple", "2" = "orange")

# Save the heatmap as a PNG file
png("correlation_heatmap_oxphos_k2_row_col.png", width = 8, height = 10, units = "in", res = 300)

# Generate the heatmap and store it in a variable
p1 <- pheatmap(
  cor_matrix,
  cluster_rows = hclust_rows,  # Use precomputed row clustering
  cluster_cols = hclust_cols,  # Use precomputed column clustering
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  color = color_palette,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 4,  # Reduce row name font size
  cutree_rows = 2,   # Split rows into 2 clusters
  cutree_cols = 2,   # Split columns into 2 clusters
  main = "Correlation Heatmap of OXPHOS Genes and SNCA, CHCHD2 & CHCHD10 (k=2 row & col)"
)

dev.off()  # Close the graphics device

# Print the heatmap object to visualize in the R console
p1

# Extract correlation values for SNCA with CHCHD2, CHCHD10
snca_chchd_corr <- cor_matrix[c("SNCA"), c("CHCHD2", "CHCHD10")]

# Extract correlation values for SNCA with OXPHOS genes
oxphos_genes <- setdiff(rownames(cor_matrix), c("SNCA", "CHCHD2", "CHCHD10"))
snca_oxphos_corr <- cor_matrix["SNCA", oxphos_genes]

# Combine correlation values into data frames
snca_corr_df <- data.frame(
  Gene = c(rep("CHCHD", length(snca_chchd_corr)), rep("OXPHOS", length(snca_oxphos_corr))),
  Correlation = c(snca_chchd_corr, snca_oxphos_corr)
)

# Perform Wilcoxon rank-sum test to compare correlation distributions
wilcox_test <- wilcox.test(
  snca_corr_df$Correlation[snca_corr_df$Gene == "CHCHD"],
  snca_corr_df$Correlation[snca_corr_df$Gene == "OXPHOS"],
  alternative = "greater"  # Test if CHCHD correlations are greater
)

# Extract p-value from the test
p_value <- wilcox_test$p.value
p_value_text <- paste("Wilcoxon p-value:", signif(p_value, 3))

# Generate the plot with p-value annotation
g1 <- ggplot(snca_corr_df, aes(x = Gene, y = Correlation, fill = Gene)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Comparison of SNCA Correlations",
       y = "Pearson Correlation with SNCA",
       x = "Gene Group") +
  scale_fill_manual(values = c("CHCHD" = "green", "OXPHOS" = "grey")) +
  annotate("text", x = 1.5, y = 0.9, 
           label = p_value_text, size = 5, color = "black") + ylim(0,1)

# Save the plot to a file
ggsave("SNCA_corr_comparison.png", plot = g1, width = 8, height = 6, dpi = 300)

# Display the plot in R
print(g1)




## Compare correlations
# Ensure oxphos_genes are present in cor_cor_matrix
available_oxphos_genes <- intersect(rownames(cor_cor_matrix), oxphos_genes)

# Warn if any OXPHOS genes are missing
missing_oxphos_genes <- setdiff(oxphos_genes, available_oxphos_genes)
if (length(missing_oxphos_genes) > 0) {
  message("Warning: The following OXPHOS genes are not found in the cor-cor matrix: ", 
          paste(missing_oxphos_genes, collapse = ", "))
}

# Extract correlation values for SNCA with available OXPHOS genes
snca_oxphos_corr <- cor_cor_matrix["SNCA", available_oxphos_genes]

# Ensure genes of interest are present
available_interest_genes <- intersect(rownames(cor_cor_matrix), genes_of_interest)

# Extract correlation values for SNCA with CHCHD2, CHCHD10
snca_chchd_corr <- cor_cor_matrix["SNCA", intersect(available_interest_genes, c("CHCHD2", "CHCHD10"))]

# Extract correlation values for SNCA with all other genes excluding interest and oxphos genes
other_genes <- setdiff(rownames(cor_cor_matrix), c(available_interest_genes, available_oxphos_genes))
snca_else_corr <- cor_cor_matrix["SNCA", other_genes]

# Prepare the data frame for plotting
snca_corr_df <- data.frame(
  GeneGroup = c(rep("CHCHD2/10", length(snca_chchd_corr)), 
                rep("OXPHOS", length(snca_oxphos_corr)), 
                rep("Else", length(snca_else_corr))),
  Correlation = c(snca_chchd_corr, snca_oxphos_corr, snca_else_corr)
)

# Perform Wilcoxon rank-sum tests
wilcox_chchd_vs_oxphos <- wilcox.test(
  snca_corr_df$Correlation[snca_corr_df$GeneGroup == "CHCHD2/10"],
  snca_corr_df$Correlation[snca_corr_df$GeneGroup == "OXPHOS"],
  alternative = "greater"
)

wilcox_chchd_vs_else <- wilcox.test(
  snca_corr_df$Correlation[snca_corr_df$GeneGroup == "CHCHD2/10"],
  snca_corr_df$Correlation[snca_corr_df$GeneGroup == "Else"],
  alternative = "greater"
)

wilcox_oxphos_vs_else <- wilcox.test(
  snca_corr_df$Correlation[snca_corr_df$GeneGroup == "OXPHOS"],
  snca_corr_df$Correlation[snca_corr_df$GeneGroup == "Else"],
  alternative = "greater"
)

# Extract p-values
p_value_oxphos <- signif(wilcox_chchd_vs_oxphos$p.value, 3)
p_value_else <- signif(wilcox_chchd_vs_else$p.value, 3)
p_value_oxphos_vs_else <- signif(wilcox_oxphos_vs_else$p.value, 3)

# Add p-values to plot text
p_value_text <- paste(
  "CHCHD2/10 vs OXPHOS p-value:", p_value_oxphos, "\n",
  "CHCHD2/10 vs Else p-value:", p_value_else, "\n",
  "OXPHOS vs Else p-value:", p_value_oxphos_vs_else
)

# Order data by mean correlation
snca_corr_df$GeneGroup <- factor(snca_corr_df$GeneGroup, 
                                 levels = names(sort(tapply(snca_corr_df$Correlation, snca_corr_df$GeneGroup, mean), decreasing = TRUE)))

# Plot the correlation data with p-value annotation
g1 <- ggplot(snca_corr_df, aes(x = GeneGroup, y = Correlation, fill = GeneGroup)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Comparison of SNCA Correlations",
       y = "Pearson Correlation with SNCA",
       x = "Gene Group") +
  scale_fill_manual(values = c("CHCHD2/10" = "green", "OXPHOS" = "grey", "Else" = "blue")) +
  annotate("text", x = 2, y = max(snca_corr_df$Correlation) * 1.1, 
           label = p_value_text, size = 4, color = "black") + ylim(0,1.2)

# Save the plot
ggsave("SNCA_corcor_comparison.png", plot = g1, width = 8, height = 6, dpi = 300)

# Display the plot in R console
print(g1)




#------------ Reviewer Q's - comparing CHCHD2 and CHCHD10 ----------

exp_dat <- as.matrix(gxdat_s@assays$RNA$data)
meta_dat <- as.data.frame(gxdat_s@meta.data)

# Filter metadata for relevant samples
meta_dat1 <- meta_dat[meta_dat$ROI %in% c("SNV") &
                        meta_dat$segment == "TH" &
                        meta_dat$Diagnosis %in% c("CTR"), ]
exp_dat1 <- exp_dat[, row.names(meta_dat1)]

# Calculate correlation matrix - total matrix
cor_matrix <- cor(t(exp_dat1), method = "pearson", use = "pairwise.complete.obs")

# Function to generate correlation plots
generate_corr_plot <- function(gene_of_interest, cor_matrix, oxphos_genes) {
  
  # Ensure oxphos_genes are present in cor_matrix
  available_oxphos_genes <- intersect(rownames(cor_matrix), oxphos_genes)
  
  # Warn if any OXPHOS genes are missing
  missing_oxphos_genes <- setdiff(oxphos_genes, available_oxphos_genes)
  if (length(missing_oxphos_genes) > 0) {
    message("Warning: The following OXPHOS genes are not found in the cor_matrix: ", 
            paste(missing_oxphos_genes, collapse = ", "))
  }
  
  # Ensure the gene of interest is present
  if (!(gene_of_interest %in% rownames(cor_matrix))) {
    stop(paste("Error: Gene", gene_of_interest, "not found in the correlation matrix."))
  }
  
  # Extract correlation values for the gene of interest with available OXPHOS genes
  interest_oxphos_corr <- cor_matrix[gene_of_interest, available_oxphos_genes]
  
  # Extract correlation values for the gene of interest with all other genes excluding itself and OXPHOS genes
  other_genes <- setdiff(rownames(cor_matrix), c(gene_of_interest, available_oxphos_genes))
  interest_else_corr <- cor_matrix[gene_of_interest, other_genes]
  
  # Prepare the data frame for plotting
  corr_df <- data.frame(
    GeneGroup = c(rep("OXPHOS", length(interest_oxphos_corr)), 
                  rep("Else", length(interest_else_corr))),
    Correlation = c(interest_oxphos_corr, interest_else_corr)
  )
  
  # Perform Wilcoxon rank-sum test
  wilcox_oxphos_vs_else <- wilcox.test(
    corr_df$Correlation[corr_df$GeneGroup == "OXPHOS"],
    corr_df$Correlation[corr_df$GeneGroup == "Else"],
    alternative = "greater"
  )
  
  # Extract p-value
  p_value_oxphos_vs_else <- signif(wilcox_oxphos_vs_else$p.value, 3)
  
  # Prepare p-value text
  p_value_text <- paste("P =", p_value_oxphos_vs_else)
  
  # Order data by mean correlation
  corr_df$GeneGroup <- factor(corr_df$GeneGroup, 
                              levels = names(sort(tapply(corr_df$Correlation, corr_df$GeneGroup, mean), decreasing = TRUE)))
  
  # Plot the correlation data with p-value annotation
  g <- ggplot(corr_df, aes(x = GeneGroup, y = Correlation, fill = GeneGroup)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.7) +
    theme_minimal() +
    labs(title = "",
         y = paste("Pearson Correlation with", gene_of_interest),
         x = "Gene Group") +
    scale_fill_manual(values = c("OXPHOS" = "grey", "Else" = "blue")) +
    annotate("text", x = 1.5, y = max(corr_df$Correlation) * 1.1, 
             label = p_value_text, size = 4, color = "black") + ylim(0,1.1)
  
  # Save the plot
  # filename <- paste0(gene_of_interest, "_correlation_comparison.png")
  # ggsave(filename, plot = g, width = 8, height = 6, dpi = 300)
  # 
  # Return the plot
  return(g)
}

# Generate plots for CHCHD2 and CHCHD10
g1 <- generate_corr_plot("CHCHD2", cor_matrix, oxphos_genes)

g2 <- generate_corr_plot("CHCHD10", cor_matrix, oxphos_genes)

# Display the plots in the R console
print(g1)
print(g2)

# arrange
arrange <- ggarrange(plotlist=list(g1,g2), nrow=2, ncol=2)
ggsave("boxplots_chchd2.10_oxphos.else_th.snv.pdf", arrange,width = 8, height = 6)




# Filter metadata for relevant samples
meta_dat1 <- meta_dat[meta_dat$ROI %in% c("SND") &
                        meta_dat$segment == "TH" &
                        meta_dat$Diagnosis %in% c("CTR"), ]
exp_dat1 <- exp_dat[, row.names(meta_dat1)]

# Calculate correlation matrix - total matrix
cor_matrix <- cor(t(exp_dat1), method = "pearson", use = "pairwise.complete.obs")

# Generate plots for CHCHD2 and CHCHD10
g1 <- generate_corr_plot("CHCHD2", cor_matrix, oxphos_genes)

g2 <- generate_corr_plot("CHCHD10", cor_matrix, oxphos_genes)

# Display the plots in the R console
print(g1)
print(g2)

# arrange
arrange <- ggarrange(plotlist=list(g1,g2), nrow=2, ncol=2)
ggsave("boxplots_chchd2.10_oxphos.else_th.snd.pdf", arrange,width = 8, height = 6)






#------------ Reviewer Q's Oxphos expression between SND and SNV ----------

# extract data matrices - counts for limma
exp_dat <- as.matrix(gxdat_s@assays$RNA$counts)
meta_dat <- as.data.frame(gxdat_s@meta.data)
colnames(meta_dat) <- make.names(colnames(meta_dat))

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
  res[[val]] <- tt
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

# Get controls
df <- df[df$contrast == "CTR",] 
df <- df[complete.cases(df),]

# Define significance categories
df$color_group <- 'grey'  # Default: all genes are grey (non-significant)

# Assign colors based on significance thresholds
df$color_group[df$adj.P.Val < 0.05 & abs(df$condtTRUE) > 0.5] <- 'red'   # Highly significant
df$color_group[df$adj.P.Val < 0.05 & abs(df$condtTRUE) <= 0.5] <- 'blue'  # p < 0.05 but FC not significant
df$color_group[df$adj.P.Val >= 0.05 & abs(df$condtTRUE) > 0.5] <- 'black' # FC significant but p not

# Set all non-OXPHOS genes to faded grey
df$color_group[!df$Gene %in% oxphos_genes] <- 'grey'

# Set labels: only for OXPHOS genes that meet significance criteria
df$label_gene <- ifelse(df$Gene %in% oxphos_genes & df$color_group %in% c("red", "blue", "black"), df$Gene, "")

# Volcano plot
plot1 <- EnhancedVolcano(df,
                         lab = df$label_gene,  # Only label significant OXPHOS genes
                         x = 'condtTRUE',  
                         y = 'adj.P.Val',  
                         pCutoff = 0.05,  
                         FCcutoff = 0.5,  
                         pointSize = 3.0,  
                         labSize = 4.0,  
                         colAlpha = ifelse(df$Gene %in% oxphos_genes, 0.75, 0.1),  # Lower alpha for non-OXPHOS genes
                         title = 'OXPHOS DEG (SNV - SND)',
                         subtitle = '',
                         caption = "Threshold: adj.P.Val < 0.05 & |FC| > 0.5",
                         legendPosition = 'right',
                         col = c('black', 'grey', 'blue', 'red'),  # Order: Black (FC only), Grey (NS), Blue (p only), Red (both)
                         drawConnectors = TRUE,  
                         widthConnectors = 1)

ggsave("volcano_oxphos_DEG.png", plot = plot1, width = 12, height = 8, dpi = 300)
