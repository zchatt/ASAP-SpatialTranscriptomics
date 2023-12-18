### public dataset processing
library(Seurat)
library(SeuratDisk)
library(data.table)

### locations ###
public_dataset_location = "/Volumes/PRJ-2022ASAPs/Bioinformics/public_datasets/" # this is the mapped project drive
seurat_obj_location = paste0(public_dataset_location,"seurat_objects/") # save location for all seurat objects
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

# v) normalise and scale data
seurat_obj_sub <- SCTransform(seurat_obj_sub)
seurat_obj_sub <- ScaleData(seurat_obj_sub)

# iv) save seurat object
save(seurat_obj, file = paste0(seurat_obj_location,"kamath_seurat.Rdata"))