# linear regression of cell-types deconvoluted (PAGE or C-SIDE) within GeoMx
library(edgeR)
library(dplyr)
library(EnhancedVolcano)

## inputs
analysis_dir <- "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis"
rdata = "/Users/zacc/USyd/spatial_transcriptomics/analysis/geomx/geomx_sep2023/analysis/geomx_sep2023_seurat.Rdata"

############################################################################################
###### Part 1: Regression analysis of C-SIDE results; brain-region x disease stage
############################################################################################
# get C-SIDE and PAGE results and format
setwd(analysis_dir)
load(rdata)
barcodes <- colnames(myRCTD@spatialRNA@counts)
norm_weights <- normalize_weights(myRCTD@results$weights)
meta <- gxdat_s@meta.data[gxdat_s@meta.data$segment.y != "TH" & gxdat_s$Diagnosis != "NTC",]
meta <- cbind(meta,norm_weights[row.names(meta),])

# setup design matrix and perform contrast
brain_regions <- c("A10","A9","A6")
cell_types <- colnames(norm_weights)
res_list <- list()
for (i in 1:length(brain_regions)){
  run_index <- meta$Brainregion %in% brain_regions[i]
  targ <- meta[run_index,]
  colnames(targ) <- make.names(colnames(targ))
  y <- t(meta[run_index,cell_types]) # C-SIDE cell-types
  #y <- y + 30 # adjust so not negative
  
  targ$DV200 <- as.numeric(targ$DV200)
  targ$Age <- as.numeric(targ$Age)
  targ$area <- as.numeric(targ$area)
  targ$IHC.score <- as.numeric(targ$IHC.score)
  targ$Dx_cont <- recode(targ$Diagnosis,"CTR" = 1,"ePD" = 2, "ILBD" = 3,"lPD" = 4 )
  
  ## Create design matrix
  design <- model.matrix(~ Dx_cont + Age + Sex + DV200 + area ,data=targ)
  colnames(design) <- make.names(colnames(design))
  v <- voom(y,design)
  vfit <- lmFit(v)
  
  # Perform LIMMA contrasts
  cont.matrix <- makeContrasts(A=Dx_cont,levels=design)
  fit2 <- contrasts.fit(vfit, cont.matrix)
  vfit2 <- eBayes(fit2)
  options(digits=3)
  
  # Select significant k-mers for each contrast and bind together
  topA <- topTable(vfit2,coef="A",number=Inf,sort.by="P")
  
  res_list[[i]] <- topA
}

# plot logFC vs p-value
pdf("volcano_limma_csidecell_full_res1.pdf")
EnhancedVolcano(toptable = res_list[[1]], x = "logFC",y = "adj.P.Val",
                lab = rownames(res_list[[1]]),pCutoff = 0.05,
                title = brain_regions[1],FCcutoff = 0,xlim = c(-0.05,0.05), ylim = c(0,10),
                subtitle = "")
dev.off()
pdf("volcano_limma_csidecell_full_res2.pdf")
EnhancedVolcano(toptable = res_list[[2]], x = "logFC",y = "adj.P.Val",
                lab = rownames(res_list[[2]]),pCutoff = 0.05,
                title = brain_regions[2],FCcutoff = 0,xlim = c(-0.05,0.05), ylim = c(0,10),
                subtitle = "")
dev.off()
pdf("volcano_limma_csidecell_full_res3.pdf")
EnhancedVolcano(toptable = res_list[[3]], x = "logFC",y = "adj.P.Val",
                lab = rownames(res_list[[3]]),pCutoff = 0.05,
                title = brain_regions[3],FCcutoff = 0,xlim = c(-0.05,0.05), ylim = c(0,10),
                subtitle = "")

dev.off()


# plot boxplots of top hits
d1 <- meta[meta$Brainregion %in% brain_regions,]

p1a <- ggplot(d1, aes(Diagnosis, y=Astro_CYP4F12, fill= Diagnosis)) +
  geom_boxplot(position=position_dodge(1)) + facet_grid( ~ Brainregion) +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif")

p1b <- ggplot(d1, aes(Diagnosis, y=SOX6_DDT, fill= Diagnosis)) +
  geom_boxplot(position=position_dodge(1)) + facet_grid( ~ Brainregion) +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif")


p1c <- ggplot(d1, aes(Diagnosis, y=Inh_OTX2_CASR, fill= Diagnosis)) +
  geom_boxplot(position=position_dodge(1)) + facet_grid( ~ Brainregion) +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif")

p1d <- ggplot(d1, aes(Diagnosis, y=CALB1_CALCR, fill= Diagnosis)) +
  geom_boxplot(position=position_dodge(1)) + facet_grid( ~ Brainregion) +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif")


p1e <- ggplot(d1, aes(Diagnosis, y=MG_TSPO_VIM, fill= Diagnosis)) +
  geom_boxplot(position=position_dodge(1)) + facet_grid( ~ Brainregion) +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif")

ggsave("limma_csidecells_a9a10a6.png", 
       ggarrange(p1a,p1b,
                 p1c,p1d,nrow=2,ncol=2), 
       device = "png")

############################################################################################
###### Part 2: Regression analysis of PAGE results; A9 v A10 for each disease stage
############################################################################################
## early changes in cell-types
meta <- gxdat_s@meta.data[gxdat_s@meta.data$segment.y != "TH" & gxdat_s$Diagnosis %in% c("CTR","ePD","ILBD","lPD"),]
meta <- cbind(meta,norm_weights[row.names(meta),])
meta$Brainregion <- factor(meta$Brainregion, levels = c("A6","A9","A10"))

brain_regions <- c("A10","A9","A6")
cell_types <- colnames(meta)[39:92]
res_list <- list()

# setup desing matrix
run_index <- meta$Brainregion %in% brain_regions
targ <- meta[run_index,]
colnames(targ) <- make.names(colnames(targ))
y <- t(meta[run_index,cell_types]) # PAGE cell-types
y <- y + 30 # adjust so not negative

targ$DV200 <- as.numeric(targ$DV200)
targ$Age <- as.numeric(targ$Age)
targ$area <- as.numeric(targ$area)
targ$IHC.score <- as.numeric(targ$IHC.score)
targ$Dx_cont <- recode(targ$Diagnosis,"CTR" = 1,"ePD" = 2, "ILBD" = 3,"lPD" = 4 )

## Create design matrix
design <- model.matrix(~ Diagnosis * Brainregion + Age + Sex + DV200 + area ,data=targ)
colnames(design) <- make.names(colnames(design))
v <- voom(y,design)
vfit <- lmFit(v)

# Perform LIMMA contrasts
cont.matrix <- makeContrasts(A=DiagnosisePD.BrainregionA9 - DiagnosisePD.BrainregionA10,
                             B=DiagnosislPD.BrainregionA9 - DiagnosislPD.BrainregionA10,
                             # C=DiagnosisILBD.BrainregionA9 - DiagnosisILBD.BrainregionA10,
                             levels=design)
fit2 <- contrasts.fit(vfit, cont.matrix)
vfit2 <- eBayes(fit2)
options(digits=3)

# Select significant k-mers for each contrast and bind together
topA <- topTable(vfit2,coef="A",number=Inf,sort.by="P")
topB <- topTable(vfit2,coef="B",number=Inf,sort.by="P")

# plot boxplots of top hits
d1 <- meta[meta$Brainregion %in% c("A9","A10"),]

p1a <- ggplot(d1, aes(Brainregion, y=CD8_Tpex.Andreatta.et.al.2021 , fill= Brainregion)) +
  geom_boxplot(position=position_dodge(1)) + facet_grid( ~ Diagnosis) +  xlab("CD8 Tpex (PAGE)")  +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif")

p1b <- ggplot(d1, aes(Brainregion, y=Treg.Andreatta.et.al.2021  , fill= Brainregion)) +
  geom_boxplot(position=position_dodge(1)) + facet_grid( ~ Diagnosis) +  xlab("CD8 Tpex (PAGE)")  +
  scale_fill_brewer(palette='viridis') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text=element_text(size=10)) + 
  stat_compare_means( label = "p.signif")

ggsave("limma_csidecells_ePDa9a10.lPDa9a10.png", 
       ggarrange(p1a,nrow=2,ncol=2), 
       device = "png")