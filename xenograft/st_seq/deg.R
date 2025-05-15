# evaluation of genes of interest and DEG in xenograft Visium data

# libs
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)


#----------- Input -----------

# working dir
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft")

# load Seurat object
seurat_obj_merge <- LoadSeuratRds( "/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft/xenograft_seurat_250325.rds")

# subset grafts
seurat_graft <- subset(seurat_obj_merge, subset = grafted_area == "TRUE")


#----------- Genes of interest expression evaluation -----------

### 1) immature neurons

# Define genes of interest
Immature_neurons <- c("TBR1","NEUROD1","DCX","STMN1")
genes_of_interest <- Immature_neurons 
genes_of_interest[!genes_of_interest %in% row.names(seurat_graft@assays$SCT$data)]

# Get SCT-normalized expression for selected genes
expr_data <- as.data.frame(as.matrix(scale(seurat_graft@assays$SCT$data)[genes_of_interest, ])) # genes x cells scaled column wise
expr_data <- t(expr_data) %>% as.data.frame()  # transpose to cells x genes
expr_data <- rownames_to_column(expr_data, var = "Cell")

# Extract metadata (Genotype and Graft), and add cell names
meta <- seurat_graft@meta.data %>%
  rownames_to_column(var = "Cell") %>%
  dplyr::select(Cell, Genotype, Graft)

# Merge expression with metadata
expr_long <- left_join(expr_data, meta, by = "Cell")

# Create group label and pivot longer for ggplot
expr_long <- expr_long %>%
  mutate(Group = paste(Genotype, Graft, sep = "_")) %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression")

# Define custom colors for groups
group_colors <- c(
  "Control_Cortex" = "lightcoral",
  "Control_VM" = "lightblue",   # light red
  "LRRK2_Cortex" = "red3",
  "PRKN_Cortex" = "pink",
  "SNCA_Cortex" = "red"
)

combined_plot <- ggplot(expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", color = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  ggtitle("Gene Expression (SCT normalized & scaled) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot

# Save
ggsave("Vln_genesinterest_immature.neurons.png", combined_plot)




### 2) mature neurons

# Define genes of interest
Mature_neurons <- c("RBFOX3","MAP2","NEFM","NEFH","SYP","DLG4")
genes_of_interest <- Mature_neurons
genes_of_interest[!genes_of_interest %in% row.names(seurat_graft@assays$SCT$data)]

# Get SCT-normalized expression for selected genes
expr_data <- as.data.frame(as.matrix(scale(seurat_graft@assays$SCT$data)[genes_of_interest, ])) # genes x cells scaled column wise
expr_data <- t(expr_data) %>% as.data.frame()  # transpose to cells x genes
expr_data <- rownames_to_column(expr_data, var = "Cell")

# Extract metadata (Genotype and Graft), and add cell names
meta <- seurat_graft@meta.data %>%
  rownames_to_column(var = "Cell") %>%
  dplyr::select(Cell, Genotype, Graft)

# Merge expression with metadata
expr_long <- left_join(expr_data, meta, by = "Cell")

# Create group label and pivot longer for ggplot
expr_long <- expr_long %>%
  mutate(Group = paste(Genotype, Graft, sep = "_")) %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression")

# Define custom colors for groups
group_colors <- c(
  "Control_Cortex" = "lightcoral",
  "Control_VM" = "lightblue",   # light red
  "LRRK2_Cortex" = "red3",
  "PRKN_Cortex" = "pink",
  "SNCA_Cortex" = "red"
)

combined_plot <- ggplot(expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", color = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  ggtitle("Gene Expression (SCT normalized & scaled) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot

# Save
ggsave("Vln_genesinterest_mature.neurons.png", combined_plot)



### 3) Dopamine neurons
Dopamine_neurons <- c("TH","SLC6A3","FOXA2","KCNJ6","NR4A2","LMX1B","SOX6","CALB1","GAD2","ALDH1A1")
genes_of_interest <- Dopamine_neurons 
genes_of_interest[!genes_of_interest %in% row.names(seurat_graft@assays$SCT$data)]

# Get SCT-normalized expression for selected genes
expr_data <- as.data.frame(as.matrix(scale(seurat_graft@assays$SCT$data)[genes_of_interest, ])) # genes x cells scaled column wise
expr_data <- t(expr_data) %>% as.data.frame()  # transpose to cells x genes
expr_data <- rownames_to_column(expr_data, var = "Cell")

# Extract metadata (Genotype and Graft), and add cell names
meta <- seurat_graft@meta.data %>%
  rownames_to_column(var = "Cell") %>%
  dplyr::select(Cell, Genotype, Graft)

# Merge expression with metadata
expr_long <- left_join(expr_data, meta, by = "Cell")

# Create group label and pivot longer for ggplot
expr_long <- expr_long %>%
  mutate(Group = paste(Genotype, Graft, sep = "_")) %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression")

# Define custom colors for groups
group_colors <- c(
  "Control_Cortex" = "lightcoral",
  "Control_VM" = "lightblue",   # light red
  "LRRK2_Cortex" = "red3",
  "PRKN_Cortex" = "pink",
  "SNCA_Cortex" = "red"
)

combined_plot <- ggplot(expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", color = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  ggtitle("Gene Expression (SCT normalized & scaled) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot

# Save
ggsave("Vln_genesinterest_DA.neurons.png", combined_plot)




### 4) PD
PD_related <- c("SNCA","PRKN","LRRK2")
genes_of_interest <- PD_related
genes_of_interest[!genes_of_interest %in% row.names(seurat_graft@assays$SCT$data)]

# Get SCT-normalized expression for selected genes
expr_data <- as.data.frame(as.matrix(scale(seurat_graft@assays$SCT$data)[genes_of_interest, ])) # genes x cells scaled column wise
expr_data <- t(expr_data) %>% as.data.frame()  # transpose to cells x genes
expr_data <- rownames_to_column(expr_data, var = "Cell")

# Extract metadata (Genotype and Graft), and add cell names
meta <- seurat_graft@meta.data %>%
  rownames_to_column(var = "Cell") %>%
  dplyr::select(Cell, Genotype, Graft)

# Merge expression with metadata
expr_long <- left_join(expr_data, meta, by = "Cell")

# Create group label and pivot longer for ggplot
expr_long <- expr_long %>%
  mutate(Group = paste(Genotype, Graft, sep = "_")) %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression")

# Define custom colors for groups
group_colors <- c(
  "Control_Cortex" = "lightcoral",
  "Control_VM" = "lightblue",   # light red
  "LRRK2_Cortex" = "red3",
  "PRKN_Cortex" = "pink",
  "SNCA_Cortex" = "red"
)

combined_plot <- ggplot(expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", color = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  ggtitle("Gene Expression (SCT normalized & scaled) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot

# Save
ggsave("Vln_genesinterest_PD.related.png", combined_plot)



#----------- Spatial plots -----------


plot_names <- unique(paste(c$Graft, seurat_graft$Case_ID, seurat_graft$Sample_Name, sep=" / "))
list_sr <- unique(seurat_graft$Sample_Name)
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialFeaturePlot(seurat_graft, features = "DCX", images = list_sr[i],pt.size.factor = 5,image.alpha = 0.5) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.3, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("DCX_Spatial.png", arrange)



plot_names <- unique(paste(seurat_graft$Graft, seurat_graft$Case_ID, seurat_graft$Sample_Name, sep=" / "))
list_sr <- unique(seurat_graft$Sample_Name)
plot_list <- list()
for (i in 1:length(list_sr)) {
  
  plot_list[[i]] <- SpatialFeaturePlot(seurat_graft, features = "MAP2", images = list_sr[i],pt.size.factor = 5,image.alpha = 0.5) + 
    theme(legend.position = "right",legend.title = element_blank(),legend.key.size = unit(0.3, 'cm')) + 
    ggtitle(plot_names[i]) + geom_point(shape=21, stroke=NA, fill=NA, size=0.5) +  theme(plot.title = element_text(size = 8))
}

arrange <- ggarrange(plotlist=plot_list, nrow=4, ncol=4, widths = c(2,2))
ggsave("MAP2_Spatial.png", arrange)



#----------- DEG -----------
library(edgeR)

### 1) control; cortical vs VM
# define count data for whole transcriptome DEG analysis
count_mat <- seurat_graft@assays$SCT$counts
meta_dat <- as.data.frame(seurat_graft@meta.data)
table(colnames(count_mat) == row.names(meta_dat))

# select controls
meta_dat <- meta_dat[meta_dat$Genotype == "Control",]
count_mat <- count_mat[,row.names(meta_dat)]
table(colnames(count_mat) == row.names(meta_dat))

#format covariates
meta_dat$run <- as.character(meta_dat$run)
meta_dat$nCount_Spatial <- as.numeric(meta_dat$nCount_Spatial)
meta_dat$Graft <- as.character(meta_dat$Graft)
meta_dat$Case_ID <- as.character(meta_dat$Case_ID)

# define condition
meta_dat$condt <- meta_dat$Graft == "Cortex"


# Create DGEList object
dge <- DGEList(count_mat, group = meta_dat$condt)

# Filter lowly expressed genes
keep <- filterByExpr(dge, design = model.matrix(~ condt + run + nCount_Spatial, data = meta_dat))
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Proceed with normalization and rest of your pipeline
dge <- calcNormFactors(dge)

design <- model.matrix(~ condt + run + nCount_Spatial, data = meta_dat)
vm <- voom(dge, design = design, plot = TRUE)
fit <- lmFit(vm, design = design)
fit <- eBayes(fit)
tt <- topTable(fit, n = Inf, adjust.method = "BH")

# expression and p-value cutoff
res <- tt[abs(tt$condtTRUE) > 0.5 & tt$adj.P.Val < 0.05,]



### Plot
res <- res[order(-abs(res$condtTRUE)),]
genes_of_interest <- row.names(res)[1:5]

# Get SCT-normalized expression for selected genes
expr_data <- as.data.frame(as.matrix(scale(seurat_graft@assays$SCT$data)[genes_of_interest,colnames(count_mat)])) # genes x cells scaled column wise
expr_data <- t(expr_data) %>% as.data.frame()  # transpose to cells x genes
expr_data <- rownames_to_column(expr_data, var = "Cell")

# Extract metadata (Genotype and Graft), and add cell names
meta <- meta %>%
  rownames_to_column(var = "Cell") %>%
  dplyr::select(Cell, Genotype, Graft)

# Merge expression with metadata
expr_long <- left_join(expr_data, meta, by = "Cell")

# Create group label and pivot longer for ggplot
expr_long <- expr_long %>%
  mutate(Group = paste(Genotype, Graft, sep = "_")) %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression")

# Define custom colors for groups
group_colors <- c(
  "Control_Cortex" = "lightcoral",
  "Control_VM" = "lightblue"
)

combined_plot <- ggplot(expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", color = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  ggtitle("Gene Expression (SCT normalized & scaled) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot

# Save
ggsave("Vln_deg_top5_ctrCortex.ctrVM.png", combined_plot)













### 2) Cortical SNCA vs Cortical Control

# define count data for whole transcriptome DEG analysis
count_mat <- seurat_graft@assays$SCT$counts
meta_dat <- as.data.frame(seurat_graft@meta.data)
table(colnames(count_mat) == row.names(meta_dat))

# select controls
meta_dat <- meta_dat[meta_dat$Genotype %in% c("Control","SNCA")  & meta_dat$Graft == "Cortex" ,]
count_mat <- count_mat[,row.names(meta_dat)]
table(colnames(count_mat) == row.names(meta_dat))

#format covariates
meta_dat$run <- as.character(meta_dat$run)
meta_dat$nCount_Spatial <- as.numeric(meta_dat$nCount_Spatial)
meta_dat$Graft <- as.character(meta_dat$Graft)
meta_dat$Case_ID <- as.character(meta_dat$Case_ID)

# define condition
meta_dat$condt <- meta_dat$Genotype == "Control"

# Create DGEList object
dge <- DGEList(count_mat, group = meta_dat$condt)

# Filter lowly expressed genes
keep <- filterByExpr(dge, design = model.matrix(~ condt + run + nCount_Spatial, data = meta_dat))
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Proceed with normalization and rest of your pipeline
dge <- calcNormFactors(dge)

design <- model.matrix(~ condt + run + nCount_Spatial, data = meta_dat)
vm <- voom(dge, design = design, plot = TRUE)
fit <- lmFit(vm, design = design)
fit <- eBayes(fit)
tt <- topTable(fit, n = Inf, adjust.method = "BH")

# expression and p-value cutoff
# res <- tt[abs(tt$condtTRUE) > 0.5 & tt$adj.P.Val < 0.05,]
res <- tt[ tt$adj.P.Val < 0.05,]


### Plot
res <- res[order(-abs(res$condtTRUE)),]

res <- res[order(res$condtTRUE),]

genes_of_interest <- row.names(res)[1:5]

# Get SCT-normalized expression for selected genes
expr_data <- as.data.frame(as.matrix(scale(seurat_graft@assays$SCT$data)[genes_of_interest,colnames(count_mat)])) # genes x cells scaled column wise
expr_data <- t(expr_data) %>% as.data.frame()  # transpose to cells x genes
expr_data <- rownames_to_column(expr_data, var = "Cell")

# Extract metadata (Genotype and Graft), and add cell names
meta <- meta %>%
  rownames_to_column(var = "Cell") %>%
  dplyr::select(Cell, Genotype, Graft)

# Merge expression with metadata
expr_long <- left_join(expr_data, meta, by = "Cell")

# Create group label and pivot longer for ggplot
expr_long <- expr_long %>%
  mutate(Group = paste(Genotype, Graft, sep = "_")) %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression")

# Define custom colors for groups
group_colors <- c(
  "Control_Cortex" = "lightcoral",
  "Control_VM" = "lightblue"
)

combined_plot <- ggplot(expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", color = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  ggtitle("Gene Expression (SCT normalized & scaled) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot

# Save
ggsave("Vln_deg_top5_ctrCortex.ctrVM.png", combined_plot)
















# 3) Dopamine neurons

# Define genes of interest
genes_of_interest <- c("WNT5A","SEMA3A","NTN1","KCNJ6","DCC",
                       "TBR1")

Immature_neurons <- c("TBR1","NEUROD1","DCX","STMN1")

Mature_neurons <- c("RBFOX3","MAP2","NEFM","NEFH","SYP","PSD95")

Dopamine_neurons <- c("TH","DAT","FOXA2","KCNJ6","NR4A2","LMX1B")

DA_neuron_subtypes <- c("SOX6","CALB1","GAD2","ALDH1A1")

genes_of_interest <- Immature_neurons 
genes_of_interest[!genes_of_interest %in% row.names(seurat_graft@assays$SCT$data)]

# Get SCT-normalized expression for selected genes
expr_data <- as.data.frame(as.matrix(scale(seurat_graft@assays$SCT$data)[genes_of_interest, ])) # genes x cells scaled column wise
expr_data <- t(expr_data) %>% as.data.frame()  # transpose to cells x genes
expr_data <- rownames_to_column(expr_data, var = "Cell")

# Extract metadata (Genotype and Graft), and add cell names
meta <- seurat_graft@meta.data %>%
  rownames_to_column(var = "Cell") %>%
  dplyr::select(Cell, Genotype, Graft)

# Merge expression with metadata
expr_long <- left_join(expr_data, meta, by = "Cell")

# Create group label and pivot longer for ggplot
expr_long <- expr_long %>%
  mutate(Group = paste(Genotype, Graft, sep = "_")) %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression")

# Define custom colors for groups
group_colors <- c(
  "Control_Cortex" = "lightcoral",
  "Control_VM" = "lightblue",   # light red
  "LRRK2_Cortex" = "red3",
  "PRKN_Cortex" = "pink",
  "SNCA_Cortex" = "red"
)

combined_plot <- ggplot(expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", color = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  ggtitle("Gene Expression (SCT normalized & scaled) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot

# Save
ggsave("Vln_genesinterest_immature.neurons.png", combined_plot)










# 2) mature neurons

# Define genes of interest
genes_of_interest <- c("WNT5A","SEMA3A","NTN1","KCNJ6","DCC",
                       "TBR1")

Immature_neurons <- c("TBR1","NEUROD1","DCX","STMN1")

Mature_neurons <- c("RBFOX3","MAP2","NEFM","NEFH","SYP","PSD95")

Dopamine_neurons <- c("TH","DAT","FOXA2","KCNJ6","NR4A2","LMX1B")

DA_neuron_subtypes <- c("SOX6","CALB1","GAD2","ALDH1A1")

genes_of_interest <- Immature_neurons 
genes_of_interest[!genes_of_interest %in% row.names(seurat_graft@assays$SCT$data)]

# Get SCT-normalized expression for selected genes
expr_data <- as.data.frame(as.matrix(scale(seurat_graft@assays$SCT$data)[genes_of_interest, ])) # genes x cells scaled column wise
expr_data <- t(expr_data) %>% as.data.frame()  # transpose to cells x genes
expr_data <- rownames_to_column(expr_data, var = "Cell")

# Extract metadata (Genotype and Graft), and add cell names
meta <- seurat_graft@meta.data %>%
  rownames_to_column(var = "Cell") %>%
  dplyr::select(Cell, Genotype, Graft)

# Merge expression with metadata
expr_long <- left_join(expr_data, meta, by = "Cell")

# Create group label and pivot longer for ggplot
expr_long <- expr_long %>%
  mutate(Group = paste(Genotype, Graft, sep = "_")) %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression")

# Define custom colors for groups
group_colors <- c(
  "Control_Cortex" = "lightcoral",
  "Control_VM" = "lightblue",   # light red
  "LRRK2_Cortex" = "red3",
  "PRKN_Cortex" = "pink",
  "SNCA_Cortex" = "red"
)

combined_plot <- ggplot(expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", color = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  ggtitle("Gene Expression (SCT normalized & scaled) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot

# Save
ggsave("Vln_genesinterest_immature.neurons.png", combined_plot)











# 1) evaluating VM and cortical genes of interest
vm_genes_interest <- c("TH","KCNJ6","SOX6","CALB1")








high_snca_a53t <- c("PPARGC1A", # https://pubmed.ncbi.nlm.nih.gov/30753527/
                    )


low_snca_a53t <- c("FABP7", # https://pubmed.ncbi.nlm.nih.gov/30753527/
                   )


# Gene Ontology Biological Process terms high in SNCA A53T
"organonitrogen compound metabolic process"  # https://pubmed.ncbi.nlm.nih.gov/30753527/


genes_interest <- c("WNT5A","TH","SEMA3A","NTN1","KCNJ6","DCC")

invitro_neuron_markers <- c("BRN2","CTIP2","TBR1","MAP2")

vm_genes_interest <- c("TH","KCNJ6","SOX6","CALB1")


## SNCA
tmp_expr <- as.data.frame(as.matrix(seurat_graft@assays$SCT$data["SNCA", ]))
colnames(tmp_expr) <- "SNCA"

# Add metadata
tmp_expr$Genotype <- seurat_graft$Genotype
tmp_expr$Graft <- seurat_graft$Graft

# Create Genotype_Graft combo
tmp_expr$Group <- paste(tmp_expr$Genotype, tmp_expr$Graft, sep = "_")

# Violin plot
plot1 <- ggplot(tmp_expr, aes(x = Group, y = SNCA)) +
  geom_violin(fill = "skyblue", scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal() +
  ggtitle("SNCA Expression (SCT normalized) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## PRKN
tmp_expr <- as.data.frame(as.matrix(seurat_graft@assays$SCT$data["PRKN", ]))
colnames(tmp_expr) <- "PRKN"

# Add metadata
tmp_expr$Genotype <- seurat_graft$Genotype
tmp_expr$Graft <- seurat_graft$Graft

# Create Genotype_Graft combo
tmp_expr$Group <- paste(tmp_expr$Genotype, tmp_expr$Graft, sep = "_")

# Violin plot
plot2 <- ggplot(tmp_expr, aes(x = Group, y = PRKN)) +
  geom_violin(fill = "skyblue", scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal() +
  ggtitle("PRKN Expression (SCT normalized) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## LRRK2
tmp_expr <- as.data.frame(as.matrix(seurat_graft@assays$SCT$data["LRRK2", ]))
colnames(tmp_expr) <- "LRRK2"

# Add metadata
tmp_expr$Genotype <- seurat_graft$Genotype
tmp_expr$Graft <- seurat_graft$Graft

# Create Genotype_Graft combo
tmp_expr$Group <- paste(tmp_expr$Genotype, tmp_expr$Graft, sep = "_")

# Violin plot
plot3 <- ggplot(tmp_expr, aes(x = Group, y = LRRK2)) +
  geom_violin(fill = "skyblue", scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal() +
  ggtitle("LRRK2 Expression (SCT normalized) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# Step 1: Define genes of interest
genes_of_interest <- c("SNCA", "PRKN", "LRRK2")

genes_of_interest <- c("PPARGC1A","FABP7","HOMER1")

# Step 2: Get SCT-normalized expression for selected genes
expr_data <- as.data.frame(as.matrix(seurat_graft@assays$SCT$data[genes_of_interest, ])) # genes x cells
expr_data <- t(expr_data) %>% as.data.frame()  # transpose to cells x genes
expr_data <- rownames_to_column(expr_data, var = "Cell")

# Step 3: Extract metadata (Genotype and Graft), and add cell names
meta <- seurat_graft@meta.data %>%
  rownames_to_column(var = "Cell") %>%
  dplyr::select(Cell, Genotype, Graft)

# Step 4: Merge expression with metadata
expr_long <- left_join(expr_data, meta, by = "Cell")

# Step 5: Create group label and pivot longer for ggplot
expr_long <- expr_long %>%
  mutate(Group = paste(Genotype, Graft, sep = "_")) %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "Gene", values_to = "Expression")

# Step 6: Plot
combined_plot <- ggplot(expr_long, aes(x = Group, y = Expression)) +
  geom_violin(fill = "skyblue", scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_minimal() +
  ggtitle("Gene Expression (SCT normalized) by Genotype x Graft") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_plot



