# scripts for running RCTD and CSIDE method to identify ctDEG

library(spacexr)
library(Seurat)
library(Matrix)
library(doParallel)
library(ggplot2)
library(plyr)
library(data.table)
library(msigdbr)
library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(EnhancedVolcano)
source("/Users/zacc/github_repo/spacexr/R/CSIDE_plots.R")
register(SerialParam())

############################################################################################
#### Inputs
############################################################################################

analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis"
rdata = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis/geomx_sep2023_seurat.Rdata"
results_folder = '/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis/RCTD_results'

# single-cell
cell_ranger_data = "/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger" # single-cell cell ranger output
filenames <- list.files("/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger", pattern="*UMAP.tsv", full.names=TRUE) # read in metadata / umap of cell IDs
cell_id <- do.call(rbind, lapply(filenames,fread))
#cell_id <- read.table(file = '/Users/zacc/USyd/spatial_transcriptomics/analysis/Kamath2023_cellranger/da_UMAP.tsv', sep = '\t', header = TRUE) # read in metadata / umap of cell IDs
cell_id <- cell_id[!cell_id$NAME == "TYPE",]

############################################################################################
###### Part 1: RCTD & CSIDE analysis
############################################################################################
# setwd
setwd(results_folder)

### Data Preprocessing and running RCTD
if(!file.exists(file.path(results_folder,'myRCTDde_gx.rds'))) {
  
  ## spatial data
  # load geomx normalised data
  load(rdata)
  # select samples
  #keep_index <- gxdat_s$Diagnosis == "CTR"
  keep_index <- gxdat_s$Diagnosis != "NTC"
  # load in counts matrix
  counts <- gxdat_s@assays$GeoMx@counts[,keep_index]  # load in counts matrix
  # create metadata matrix 
  meta <- gxdat_s@meta.data[keep_index,]
  # confirm names
  table(colnames(counts) == row.names(meta))
  # load in coordinates
  coords = as.data.frame(gxdat_s@reductions$umap@cell.embeddings[keep_index,]) # we are using UMAP coordinates as dummy variables until we obtain DSP ROI locations
  colnames(coords) <- c("imagerow","imagecol")
  # process counts
  nUMI <- colSums(counts) # In this case, total counts per spot
  puck <- SpatialRNA(coords, counts, nUMI)
  barcodes <- colnames(puck@counts) # pucks to be used (a list of barcode names). 
  plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')
  
  # #single-cell data
  sc_obj <- Read10X(cell_ranger_data) # single-cell cell ranger output
  # subset for cell-types
  sc_obj <- sc_obj[,cell_id$NAME]
  cell_types <- setNames(cell_id$Cell_Type, cell_id$NAME)
  #cell_types <- as.factor(cell_id$Cell_Type) # convert to factor data type
  cell_types <- as.factor(cell_types)
  nUMI <- colSums(sc_obj)
  
  ## create reference and run RCTD
  reference <- Reference(sc_obj, cell_types, nUMI)
  myRCTD <- create.RCTD(puck, reference, max_cores = 8)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  saveRDS(myRCTD,file.path(results_folder,'myRCTDde_gx.rds'))
} else {
  myRCTD <- readRDS(file.path(results_folder,'myRCTDde_gx.rds'))
}

### Exploring the full mode results
barcodes <- colnames(myRCTD@spatialRNA@counts)
norm_weights <- normalize_weights(myRCTD@results$weights)
myRCTD@config$doublet_mode <- "full"
print(head(norm_weights)) 

p1 <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'SOX6_DDT'], 
                           title ='plot of SOX6_DDT weights', size = 0.5) 

ggsave("cside_1_sox6dtt.png", p1, device = "png")

p2 <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'CALB1_RBP4'], ylimit = c(0,0.5), 
                           title ='plot of CALB1_RBP4 weights', size = 0.5) 

ggsave("cside_2_calb1rbp4.png", p2, device = "png")


############################################################################################
###### Part X: Run C-SIDE two regions
############################################################################################
meta <- gxdat_s@meta.data[gxdat_s@meta.data$segment.y == "TH" & gxdat_s@meta.data$Diagnosis == "CTR" & gxdat_s@meta.data$Brainregion %in% c("A9","A10"),]
barcodes_run <- row.names(meta)
myRCTD@config$max_cores <- 8
cell_type_threshold = 0.001
agg_cell_types <- aggregate_cell_types(myRCTD, barcodes_run,doublet_mode = F)
cell_run <- names(agg_cell_types)[agg_cell_types > cell_type_threshold]

explanatory.variable <- as.numeric(as.factor(meta$Brainregion)) - 1
names(explanatory.variable) <-  barcodes_run 

myRCTD_ctr.a9a10 <- run.CSIDE.single(myRCTD,
                           explanatory.variable,
                           cell_types = cell_run,
                           cell_type_threshold = cell_type_threshold,
                           fdr = 0.25, 
                           doublet_mode = FALSE) 

save(myRCTD_ctr.a9a10,file="myRCTD_ctr.a9a10")

## CSIDE results
all_gene_list <- myRCTD_ctr.a9a10@de_results$all_gene_list
sig_gene_list <- lapply(all_gene_list, function(x) x[x$p_val < 0.05,])
sig_gene_list <- lapply(sig_gene_list, function(x) x[order(x$p_val),])
lapply(sig_gene_list,head)

## Plots
# barplot for n sig genes x cell-type
x = unlist(lapply(sig_gene_list, nrow))

dplot <- as.data.frame(cbind(cell_type=names(x),sig_genes=x))
dplot$sig_genes <- as.numeric(dplot$sig_genes) 

g2 <- ggplot(dplot, aes(x=reorder(cell_type,-x),y=x)) +
  geom_bar(stat="identity", width=0.7)+  
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("Cell Type") + ylab("Sig CtDEG's (n)") 

# save plots
ggsave("ctDEG_results_bar_th.ctr.a10a9rn.png", 
       ggarrange(g2,nrow=2,ncol=1)
       , device = "png")


# volcano plots
dat_tmp <- as.data.frame(all_gene_list$SOX6_DDT)
v1 <- EnhancedVolcano(toptable = dat_tmp , x = "log_fc",y = "p_val",
                lab = rownames(dat_tmp),pCutoff = 0.05, FCcutoff = 0.5,
                title = "volcano_ctDEG_geomx_th.ctr.a10a9rn",subtitle = "")

ggsave("ctDEG_results_volcano_th.ctr.a10a9rn.png", 
       ggarrange(v1,nrow=1,ncol=1)
       , device = "png")

# define genesets
all_gene_sets = msigdbr(species = "human")
# using "H: hallmark gene sets"
# using "C2: curated gene sets"
# using "C5: ontology gene sets"
# using "C8: cell type signature gene sets"
msigdbr_df <- all_gene_sets %>%
  dplyr::filter(gs_cat == "C5")

# create list for fgseas
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# select res from list, perform FGSEA and plot results
#for (i in 1:length(all_gene_list)){
  i = which(names(sig_gene_list) == "SOX6_DDT" )
  print(i)
  res <- as.data.frame(sig_gene_list[[i]])
  name1 <- names(sig_gene_list)[i]
  gene_stats <- res$log_fc
  names(gene_stats) <- row.names(res)
  
  # FGSEA
  fgseaRes <- fgsea(pathways = msigdbr_list,
                    stats    = gene_stats,
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes <- fgseaRes[order(pval), ]
  head(fgseaRes)
  
  # plot
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  pdf(paste0("FGSEA_C5celltype_cside_th.ctr.a10a9rn",name1,".pdf"))
  print(plotGseaTable(msigdbr_list[topPathways], gene_stats, fgseaRes,
                      gseaParam=0.5, pathwayLabelStyle=list(size=6, color="black")))
  dev.off()
#}

############################################################################################
###### Part X: Run C-SIDE multiple regions
############################################################################################
# ## Running CSIDE - regions/ multiple contrasts
meta <- gxdat_s@meta.data[gxdat_s@meta.data$segment.y != "TH" & gxdat_s@meta.data$Diagnosis == "CTR",]
barcodes_run <- row.names(meta)
cell_type_threshold = 1
agg_cell_types <- colSums(myRCTD@results$weights[barcodes_run, ])
cell_run <- names(agg_cell_types)[agg_cell_types > (cell_type_threshold)]
myRCTD@config$max_cores <- 8

A9 <- barcodes_run[which(meta$Brainregion == "A9")]
A10 <- barcodes_run[which(meta$Brainregion == "A10")]
RN <- barcodes_run[which(meta$Brainregion == "RN")]
region_list <- list(A9,A10,RN)

X <- build.designmatrix.regions(myRCTD,region_list)
barcodes2 <- rownames(X)
myRCTD_ctr.fulla9a10rn <- run.CSIDE(myRCTD, X, barcodes2, cell_run , 
                              test_mode = 'categorical', gene_threshold = 5e-5, 
                              doublet_mode = F, cell_type_threshold = cell_type_threshold, fdr = 0.25, weight_threshold = 0, 
                              log_fc_thresh = 0)

save(myRCTD_ctr.fulla9a10rn,file="myRCTD_ctr.fulla9a10rn")

## CSIDE results
all_gene_list <- myRCTD_ctr.fulla9a10rn@de_results$all_gene_list
sig_gene_list <- lapply(all_gene_list, function(x) x[x$p_val_best < 0.05,])
sig_gene_list <- lapply(sig_gene_list, function(x) x[order(x$p_val_best),])
lapply(sig_gene_list,head)

## Plots
# barplot for n sig genes x cell-type
x = unlist(lapply(sig_gene_list, nrow))

dplot <- as.data.frame(cbind(cell_type=names(x),sig_genes=x))
dplot$sig_genes <- as.numeric(dplot$sig_genes) 

g2 <- ggplot(dplot, aes(x=reorder(cell_type,-x),y=x)) +
  geom_bar(stat="identity", width=0.7)+  
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("Cell Type") + ylab("Sig CtDEG's (n)") 

# save plots
ggsave("ctDEG_results_volcano_full.ctr.a10a9rn.png", 
       ggarrange(g2,nrow=2,ncol=1)
       , device = "png")


# volcano plots
dat_tmp <- as.data.frame(all_gene_list$SOX6_DDT)

EnhancedVolcano(toptable = dat_tmp , x = "log_fc_best",y = "p_val_best",
                lab = rownames(dat_tmp),pCutoff = 0.05, FCcutoff = 0.5,
                title = "volcano_ctDEG_geomx_CTR",subtitle = "")

# define genesets
all_gene_sets = msigdbr(species = "human")
# using "H: hallmark gene sets"
# using "C2: curated gene sets"
# using "C5: ontology gene sets"
# using "C8: cell type signature gene sets"
msigdbr_df <- all_gene_sets %>%
  dplyr::filter(gs_cat == "C5")

# create list for fgseas
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# select res from list, perform FGSEA and plot results
for (i in 1:length(all_gene_list)){
  i = which(names(sig_gene_list) == "SOX6_DDT" )
  print(i)
  res <- as.data.frame(sig_gene_list[[i]])
  name1 <- names(sig_gene_list)[i]
  gene_stats <- res$log_fc
  names(gene_stats) <- row.names(res)
  
  # FGSEA
  fgseaRes <- fgsea(pathways = msigdbr_list,
                    stats    = gene_stats,
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes <- fgseaRes[order(pval), ]
  head(fgseaRes)
  
  # plot
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  pdf(paste0("FGSEA_C5celltype_cside_fullroi.ctr.a10a9rn",name1,".pdf"))
  print(plotGseaTable(msigdbr_list[topPathways], gene_stats, fgseaRes,
                      gseaParam=0.5, pathwayLabelStyle=list(size=6, color="black")))
  dev.off()
}

# quantitative expression plot
barcodes <- myRCTD_ctr.fulla9a10rn@internal_vars_de$barcodes
gene <- 'ATP6V0B'
Y <- myRCTD@spatialRNA@counts[gene, barcodes]
Yn <- Y / myRCTD_ctr.fulla9a10rn@spatialRNA@nUMI[barcodes]
dc_prop <- myRCTD_ctr.fulla9a10rn@internal_vars_de$my_beta[,'SOX6_DDT']
thresh <- median(dc_prop)
reg <- myRCTD_ctr.fulla9a10rn@internal_vars_de$X2[barcodes,2]
pred <- predict_CSIDE_all(myRCTD_ctr.fulla9a10rn, gene)


gene_fits <- myRCTD_ctr.fulla9a10rn@de_results$gene_fits
cell_types <- "SOX6_DDT"
my_beta <- myRCTD_ctr.fulla9a10rn@internal_vars_de$my_beta
cell_type_ind <- which(myRCTD_ctr.fulla9a10rn@internal_vars_de$cell_types == "SOX6_DDT")
X2_mat <- myRCTD_ctr.fulla9a10rn@internal_vars_de$X2[rownames(my_beta),]

pred_ct <- predict_CSIDE(cell_type_ind, gene_fits, gene, X2_mat)



ct_thresh <- 0
myRCTDde <- myRCTD_ctr.fulla9a10rn
cur_cell_types <- "SOX6_DDT"
p_df <- get_quant_df(myRCTDde, myRCTDde@de_results$gene_fits, myRCTDde@internal_vars_de$cell_types, 
                     cur_cell_types, gene, multi_region = T, prop_thresh = ct_thresh)

table(row.names(p_df) == row.names(meta[barcodes,]))

d1 <- cbind(meta[barcodes,],gene_exp = p_df$pred)
d1$Brainregion <- factor(d1$Brainregion, levels = c("A9","A10","RN"))
p1a <- ggplot(d1, aes(x=Brainregion, y=gene_exp, fill=Diagnosis)) +
  geom_boxplot(position=position_dodge(1),facet.by = "Brainregion", short.panel.labs = FALSE) + 
  scale_fill_brewer(palette='viridis') + theme_classic() + ylab(paste0(gene," expression")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) +
  stat_compare_means( label = "p.signif")


d2 <- cbind(meta[barcodes,],gene_exp = Y )
d2$Brainregion <- factor(d2$Brainregion, levels = c("A9","A10","RN"))
p1b <- ggplot(d2, aes(x=Brainregion, y=gene_exp, fill=Diagnosis)) +
  geom_boxplot(position=position_dodge(1),facet.by = "Brainregion", short.panel.labs = FALSE) + 
  scale_fill_brewer(palette='viridis') + theme_classic() + ylab(paste0(gene," expression")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) +
  stat_compare_means( label = "p.signif")

# save plots
ggsave("boxplot_ATP6V0B_full.ctr.a10a9rn.png", 
       ggarrange(p1b ,nrow=2,ncol=1), device = "png")



############################################################################################
###### Part X: LIMMA on norm counts
############################################################################################
## linear regression
barcodes_run <- barcodes[gxdat_s@meta.data$segment.y == "TH" & gxdat_s@meta.data$Diagnosis == "CTR" & gxdat_s@meta.data$Brainregion %in% c("A9","A10")]
targ <- gxdat_s@meta.data[barcodes_run,]
y <- gxdat_s@assays$GeoMx@data[,barcodes_run]

# Create design matrix 
targ$DV200 <- as.numeric(targ$DV200)
targ$Age <- as.numeric(targ$Age)
targ$area <- as.numeric(targ$area)
targ$IHC.score <- as.numeric(targ$IHC.score)
#targ$Dx_cont <- recode(targ$Diagnosis,"CTR" = 1,"ePD" = 2, "ILBD" = 3,"lPD" = 4 )

## Create design matrix
design <- model.matrix(~ Brainregion + Age + Sex + DV200 + area ,data=targ)
colnames(design) <- make.names(colnames(design))
v <- voom(y,design)
vfit <- lmFit(v)

# Perform LIMMA contrasts
cont.matrix <- makeContrasts(A=BrainregionA9,levels=design)
fit2 <- contrasts.fit(vfit, cont.matrix)
vfit2 <- eBayes(fit2)
options(digits=3)

topA <- topTable(vfit2,coef="A",number=Inf,sort.by="P")

# volcano plots
dat_tmp <- as.data.frame(topA)

EnhancedVolcano(toptable = dat_tmp, y = "adj.P.Val",x = "logFC",
                lab = rownames(dat_tmp),pCutoff = 0.05, FCcutoff = 0.5,
                title = "volcano_ctDEG_geomx_CTR_TH_a10.a9",subtitle = "")


# define genesets
all_gene_sets = msigdbr(species = "human")
# using "H: hallmark gene sets"
# using "C2: curated gene sets"
# using "C5: ontology gene sets"
# using "C8: cell type signature gene sets"
msigdbr_df <- all_gene_sets %>%
  dplyr::filter(gs_cat == "C5")

# create list for fgseas
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# select res from list, perform FGSEA and plot results
for (i in 1:length(all_gene_list)){
  res <- as.data.frame(all_gene_list[[i]])
  name1 <- names(all_gene_list)[i]
  gene_stats <- res$log_fc_best
  names(gene_stats) <- row.names(res)
  
  # FGSEA
  fgseaRes <- fgsea(pathways = msigdbr_list,
                    stats    = gene_stats,
                    minSize  = 15,
                    maxSize  = 500,
                    scoreType = "pos")
  fgseaRes <- fgseaRes[order(pval), ]
  head(fgseaRes)
  
  # plot
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  pdf(paste0("FGSEA_C5celltype_cside_fullroi.ctr.a10a9rn",name1,".pdf"))
  print(plotGseaTable(msigdbr_list[topPathways], gene_stats, fgseaRes,
                      gseaParam=0.5, pathwayLabelStyle=list(size=6, color="black")))
  dev.off()
}

























class_num <- rep(0, length(barcodes)); names(class_num) <- barcodes
class_num[region_middle] <- 1; class_num[region_right] <- 2

p3  <- plot_class(myRCTD@spatialRNA, barcodes, factor(class_num),  title ='plot of regions')
# ggsave("cside_3_regionsplit.png", p3, device = "png")

#
## build design matrix
regions_list <- list()
contrast <- mapvalues(gxdat_s@meta.data$run_name,run_names,list_ID)
for (i in 1:length(unique(contrast))) {
  regions_list[[i]] <- barcodes[contrast == unique(contrast)[i]]
}
regions_list[[i+1]] <- NA
#
X <- build.designmatrix.regions(myRCTD,regions_list)
barcodes <- rownames(X)
myRCTD <- run.CSIDE(myRCTD, X, barcodes,cell_types = cell_types,
                    cell_type_threshold = 10,
                    fdr = 0.05,
                    doublet_mode = FALSE,
                    weight_threshold = 0.1,
                    test_mode = 'categorical')
#
saveRDS(myRCTD,file.path(results_folder,'myRCTDde_visium_SN_regions.rds'))
#
#CSIDE results
all_gene_list <- myRCTD@de_results$all_gene_list
sig_gene_list <- lapply(all_gene_list, function(x) x[x$p_val_best < 0.05,])

sig_gene_list <- lapply(all_gene_list, function(x) x[order(x$p_val_best),])






## Running CSIDE - single
# select just cobtrol samples
myRCTD@config$max_cores <- 2
myRCTD@config$doublet_mode = 'full'

explanatory.variable <- as.numeric(as.factor(mapvalues(gxdat_s@meta.data$run_name,run_names,list_ID))) - 1
names(explanatory.variable) <-  barcodes

myRCTD <- run.CSIDE.single(myRCTD,
                           explanatory.variable,
                           cell_type_threshold = 10,
                           fdr = 0.05, 
                           doublet_mode = FALSE,
                           weight_threshold = 0.1) 

saveRDS(myRCTD,file.path(results_folder,'myRCTDde_gx.rds'))

make_all_de_plots(myRCTD,results_folder)

## CSIDE results
all_gene_list <- myRCTD@de_results$all_gene_list
sig_gene_list <- lapply(all_gene_list, function(x) x[x$p_val < 0.05,])
sig_gene_list <- lapply(all_gene_list, function(x) x[order(x$p_val),])

