### public dataset processing
library(Seurat)
library(SeuratDisk)
library(data.table)
library(SpatialExperiment)
library(SingleCellExperiment)
library(WeberDivechaLCdata)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(anndata)
library(harmony)

### locations ###
public_dataset_location = "/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/" # this is the mapped project drive
setwd(public_dataset_location)

##################################################################################################################
### 1. Kamath et al single-cell data from SNpc - https://www.nature.com/articles/s41593-022-01061-1
##################################################################################################################
cell_ranger_data = paste0(public_dataset_location,"kamath_2023/") # single-cell cell ranger output

# i) load single-cell data
sc_obj <- Read10X(cell_ranger_data) # single-cell cell ranger output

# ii) create seurat object
seurat_obj <- CreateSeuratObject(sc_obj)
rm(sc_obj)

# iii) update metadata
tsv <- read.delim(paste0(cell_ranger_data,"METADATA_PD.tsv")) # metadata
tsv <- tsv[-1,]
row.names(tsv) <- tsv$NAME

seurat_obj <- AddMetaData(object = seurat_obj, metadata = tsv)

umap_file <- list.files(cell_ranger_data, pattern="*UMAP.tsv", full.names=TRUE) # read in umap of cell IDs
cell_id <- do.call(rbind, lapply(umap_file,fread))
cell_id <- cell_id[!cell_id$NAME == "TYPE",]
row.names(cell_id) <- cell_id$NAME

seurat_obj <- AddMetaData(object = seurat_obj, metadata = cell_id)

# iv) select control subjects
seurat_obj_sub <- subset(seurat_obj, subset = Status == 'Ctrl'  )
seurat_obj_sub <- subset(seurat_obj_sub, subset = organ__ontology_label == "substantia nigra pars compacta")
cell_types <- names(table(seurat_obj_sub@meta.data$Cell_Type))
seurat_obj_sub <- subset(seurat_obj_sub, subset = Cell_Type %in% cell_types )

# v) select proportion of cells
cell_types = unique(seurat_obj_sub@meta.data$Cell_Type)

res <- list()
for (i in 1:length(cell_types)){
  tmp <- seurat_obj_sub@meta.data[seurat_obj_sub@meta.data$Cell_Type == cell_types[i],]
  res[[i]] <- row.names(tmp)[sample(nrow(tmp), round(nrow(tmp) * 0.2,0))]
}

seurat_obj_sub2 <- seurat_obj_sub[,unlist(res)]

# vi) normalise and scale data
seurat_obj_sub2 <- SCTransform(seurat_obj_sub2)
seurat_obj_kamath_sub <- ScaleData(seurat_obj_sub2)

# vii) save seurat object
save(seurat_obj_kamath_sub, file = paste0(public_dataset_location,"kamath_seurat_sub.Rdata"))

##################################################################################################################
### 2. Siletti et al single-cell data from SNpc, re-processed by Carmen Abaurre 
##################################################################################################################
# i) read in hd4
tmp <- read_h5ad("/Users/zacc/USyd/spatial_transcriptomics/data/public_datasets/for_Zac/integration-kamath-siletti.h5ad", backed = NULL)

# ii) extract counts and metadata
exp_dat <- data.frame(t(as.matrix(tmp$X)))
cell_meta <- as.data.frame(tmp$obs)

# iii) create seurat object
seurat_silett_re <- CreateSeuratObject(exp_dat, project = "Siletti_recode", assay = "RNA", meta.data = cell_meta)

# iv) normalise and scale data
seurat_silett_re <- SCTransform(seurat_silett_re)
seurat_silett_re <- ScaleData(seurat_silett_re)

# save
save(seurat_silett_re, file = paste0(public_dataset_location,"siletti_re_seurat.Rdata"))


##################################################################################################################
### 3. Webber et al single-cell data from LC
##################################################################################################################
### LC snRNAseq data from Weber et al. ###
# i) get LC datasets
sce <- WeberDivechaLCdata_singleNucleus()

# ii) convert Ensemble ID to gene name
gene_name = mapIds(org.Hs.eg.db,
                   keys=row.names(sce), 
                   column="SYMBOL",
                   keytype="ENSEMBL",
                   multiVals="first")

sce_sub <- sce[!is.na(gene_name),]
row.names(sce_sub) <- gene_name[!is.na(gene_name)]

# iii) extract counts and metadata
exp_dat <- data.frame(as.matrix(logcounts(sce_sub)))
cell_meta <- as.data.frame(colData(sce_sub))

# iii) create seurat object
seurat_webber <- CreateSeuratObject(exp_dat, project = "Siletti_recode", assay = "RNA", meta.data = cell_meta)

# iv) normalise and scale data
seurat_webber <- SCTransform(seurat_webber)
seurat_webber <- ScaleData(seurat_webber)

# save
save(seurat_webber, file = paste0(public_dataset_location,"webber_re_seurat.Rdata"))



##################################################################################################################
### 4. Create a merged seurat object with Siletti, Kamath and Webber snRNAseq datasets
##################################################################################################################
# Performing integration on datasets normalized with SCTransform
# i) combine seurat lists
seurat_list <- c(seurat_obj_kamath_sub,seurat_silett_re,seurat_webber)

# ii) get variable features
seurat_list <- lapply(X = seurat_list, function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)) 

# iii) select features that are repeatedly variable across datasets for integration 
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

# iv) find anchors and intergrate data
merge.anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = 'SCT', anchor.features = features)
merge.combined.sct <- IntegrateData(anchorset = merge.anchors, normalization.method = 'SCT')

# v) Run the standard workflow for visualization and clustering
merge.combined.sct <- ScaleData(merge.combined.sct, verbose = FALSE)
merge.combined.sct <- RunPCA(merge.combined.sct, npcs = 30, verbose = FALSE)
merge.combined.sct <- RunUMAP(merge.combined.sct, reduction = "pca", dims = 1:30)
merge.combined.sct <- RunHarmony(merge.combined.sct, group.by.vars = c("dataset_merge") )
merge.combined.sct <- FindNeighbors(merge.combined.sct, reduction = "harmony", dims = 1:30)
merge.combined.sct <- FindClusters(merge.combined.sct, resolution = 0.5)

# vi) add metadata 
merge.combined.sct@meta.data$dataset_merge <- c(rep("Kamath",36930),rep("Siletti",2031),rep("Webber",20191))
merge.combined.sct@meta.data$dataset_merge <- as.factor(merge.combined.sct@meta.data$dataset_merge)
merge.combined.sct@meta.data$cell_type_merge <- c(merge.combined.sct@meta.data$Cell_Type[!is.na(merge.combined.sct@meta.data$Cell_Type)],
                                                  as.character(cell_meta$label_merged))


# vii) Visualization
pa <- DimPlot(merge.combined.sct, reduction = "umap")
ggsave("umap_merge.kam.sil.web.humancells.png", pa)
p1 <- DimPlot(merge.combined.sct, reduction = "umap", group.by = "dataset_merge") 
ggsave("umap_merge.kam.sil.web.humancells_dataset.png", p1)
p2 <- DimPlot(merge.combined.sct, reduction = "umap", split.by = "dataset_merge")
ggsave("umap_merge.kam.sil.web.humancells_dataset_split.png", p2)
p3 <- DimPlot(merge.combined.sct, reduction = "umap", group.by = 'cell_type_merge') +  guides(color = guide_legend(override.aes = list(size=4), ncol=3) )
ggsave("umap_merge.kam.sil.web.humancells_cell_labels.png", p3)
p3b <- DimPlot(merge.combined.sct, reduction = "umap", group.by = 'cell_type_merge') + ggplot2::theme(legend.position = "none")
ggsave("umap_merge.kam.sil.web.humancells_cell.png", p3b)


NE <- which(merge.combined.sct@meta.data$cell_type_merge == c("NE"))
p_ne<- DimPlot(merge.combined.sct, label=F, group.by="cell_type_merge", 
        cells.highlight= list(NE), cols.highlight = c("red"), cols= "grey")  + ggtitle("NE")
ggsave("umap_merge.kam.sil.web.humancells_NE.png", p_ne)


HT <- which(merge.combined.sct@meta.data$cell_type_merge == c("5HT"))
p_HT<- DimPlot(merge.combined.sct, label=F, group.by="cell_type_merge", 
               cells.highlight= list(HT), cols.highlight = c("red"), cols= "grey") + ggtitle("5-HT")
ggsave("umap_merge.kam.sil.web.humancells_5HT.png",p_HT)

toMatch <- c("5HT", "CALB1", "SOX6","Ex","NE","Inh","GAD2","inhib","neuron","excit")
neuron <- grep(paste(toMatch,collapse="|"), merge.combined.sct@meta.data$cell_type_merge)
p_neuron <- DimPlot(merge.combined.sct, label=F, group.by="cell_type_merge", 
               cells.highlight= list(neuron), cols.highlight = c("red"), cols= "grey")  + ggtitle("Neuron")
ggsave("umap_merge.kam.sil.web.humancells_neuron.png", p_neuron)

toMatch <- c("CALB1", "SOX6","GAD2")
DA <- grep(paste(toMatch,collapse="|"), merge.combined.sct@meta.data$cell_type_merge)
p_DA <- DimPlot(merge.combined.sct, label=F, group.by="cell_type_merge", 
                    cells.highlight= list(DA), cols.highlight = c("red"), cols= "grey")  + ggtitle("DA")
ggsave("umap_merge.kam.sil.web.humancells_DA.png", p_DA)

toMatch <- c("GAD2")
GAD2 <- grep(paste(toMatch,collapse="|"), merge.combined.sct@meta.data$cell_type_merge)
p_GAD2 <- DimPlot(merge.combined.sct, label=F, group.by="cell_type_merge", 
                cells.highlight= list(GAD2), cols.highlight = c("red"), cols= "grey")  + ggtitle("GAD2")
ggsave("umap_merge.kam.sil.web.humancells_GAD2.png", p_GAD2)

toMatch <- c("SOX6")
SOX6 <- grep(paste(toMatch,collapse="|"), merge.combined.sct@meta.data$cell_type_merge)
p_SOX6 <- DimPlot(merge.combined.sct, label=F, group.by="cell_type_merge", 
                cells.highlight= list(SOX6), cols.highlight = c("red"), cols= "grey")  + ggtitle("SOX6")
ggsave("umap_merge.kam.sil.web.humancells_SOX6.png", p_SOX6)


toMatch <- c("CALB1")
CALB1 <- grep(paste(toMatch,collapse="|"), merge.combined.sct@meta.data$cell_type_merge)
p_CALB1 <- DimPlot(merge.combined.sct, label=F, group.by="cell_type_merge", 
                cells.highlight= list(CALB1), cols.highlight = c("red"), cols= "grey")  + ggtitle("CALB1")
ggsave("umap_merge.kam.sil.web.humancells_CALB1.png", p_CALB1)

# counts per cell-type
VlnPlot(merge.combined.sct, features = c("GAPDH"), pt.size = 0.5,group.by = "dataset_merge",assay = "SCT") + 
  ggplot2::theme(legend.position = "none")

VlnPlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"),
        pt.size = 0.2, ncol = 4)

raw_sum <- c(colSums(merge.combined.sct[["RNA"]]$counts.1),
             colSums(merge.combined.sct[["RNA"]]$counts.2),
             colSums(merge.combined.sct[["RNA"]]$counts.3))


counts_sc <- apply(sc_obj[["integrated"]]$data,2,as.integer)

counts_sc <- sc_obj[["SCT"]]$counts
dplot <- as.data.frame(cbind(y = colSums(counts_sc),
                             x= as.character(merge.combined.sct@meta.data$dataset_merge)))
dplot$y <- as.numeric(dplot$y)
ggplot(dplot, aes(x=x, y=y)) + 
  geom_boxplot()

                             x= as.character(merge.combined.sct@meta.data$dataset_merge)))
dplot$y <- as.numeric(dplot$y)
ggplot(dplot, aes(x=x, y=y)) + 
  geom_boxplot()


sc_obj[["SCT"]]$counts

# save
save(merge.combined.sct, file = paste0(public_dataset_location,"merged_kam.sil.web_seurat.Rdata"))
