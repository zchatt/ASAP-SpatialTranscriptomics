# Validation of DEG within substantia nigra of parkinsons disease

library(readxl)

# repo location
git_repo_location = "/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/deg_markers/tables/"

#1) Ferraro et al., Cells, 2022
# Limited overlap in DEGs from published microarray studies (n=9) is caused by
# heterogeneity in the neurodegeneration and technical artefacts. 
# Meta-analysis of SNpc microarray datasets controlling for cell-type
# report LM results with a cell-type aware and unaware model
tmp <- read_excel(paste0(git_repo_location,"CorrectingForCellPropPDsn_SuppTableS3.xlsx"),sheet = 1, skip = 2 )

tmp2 <- tmp[tmp$is_cell_unaware_LMM_DEG == "yes",]
tmp2 <- tmp2[,c("Gene","cell_unaware_LMM_estimate")]
colnames(tmp2) <- c("Gene","Change")
tmp2 <- tmp2[tmp2$Gene != "NA",]
tmp2 <- as.data.frame(tmp2[complete.cases(tmp2),])
row.names(tmp2) <- tmp2$Gene
tmp2$Source <- rep("Ferraro.2022.cellunaware",nrow(tmp2))
ferraro_deg_cellunaware <- tmp2

tmp2 <- tmp[tmp$is_cell_aware_LMM_DEG == "yes" ,]
tmp2 <- tmp2[,c("Gene","cell_aware_LMM_estimate")]
colnames(tmp2) <- c("Gene","Change")
tmp2 <- tmp2[tmp2$Gene != "NA",]
tmp2 <- as.data.frame(tmp2[complete.cases(tmp2),])
row.names(tmp2) <- tmp2$Gene
tmp2$Source <- rep("Ferraro.2022.cellaware",nrow(tmp2))
ferraro_deg_cellaware <- tmp2

ferraro_deg_list <- list(ferraro_deg_cellunaware,ferraro_deg_cellaware)
names(ferraro_deg_list) <- c("Ferraro.2022.cellunaware","Ferraro.2022.cellaware")

#2) Smajić et al., Brain, 2022
# snRNAseq of PD and CTR SNpc
# ID population of DaN with high CADPS2
# increased number of activated microglia
# increased number of astrogliosis

# Differential up-regulated genes for cell-types trajectories and in PD for each cell-type
deg_names <- c(#"Microglia_TG",
               "Microglia_IPD_DEG",
               #"Astrocyte_TG",
               "Astrocyte_IPD_DEG",
               #"Oligodendrocyte_TG",
                "Oligodendrocytes_IPD_DEG",
                "Inhibitory_IPD_DEG",
                "Excitatory_IPD_DEG",
                "GABA_IPD_DEG",
                "DaNs_IPD_DEG",
                "CADPS2high_IPD_DEG",
                "OPC_IPD_DEG",
                "Ependymal_IPD_DEG",
                "Pericytes_IPD_DEG",
                "Endothelial_cell_IPD_DEG")
smajic_deg <- list()
tabs <- c(2,4,6:15)
for (k in 1:length(tabs)){
  i = tabs[k]
  tmp <- read_excel(paste0(git_repo_location,"Supplementary Table 8.xlsx"),sheet = i, skip = 3 )
  tmp <- tmp[,c(2,4)]
  colnames(tmp) <- c("Gene","Change")
  tmp <- tmp[tmp$Gene != "NA",]
  tmp <- as.data.frame(tmp[complete.cases(tmp),])
  row.names(tmp) <- tmp$Gene
  tmp$Source <- rep(paste0("Smajić.2022.", deg_names[k]),nrow(tmp))
  smajic_deg[[k]] <- tmp
  
}
names(smajic_deg) <- paste0("Smajić.2022.", deg_names)

# for trajectory, using Morans Index
deg_names <- c("Microglia_TG",
               "Astrocyte_TG",
               "Oligodendrocyte_TG")
smajic_traj <- list()
tabs <- c(1,3,5)
for (k in  1:length(tabs)){
  i = tabs[k]
  tmp <- read_excel(paste0(git_repo_location,"Supplementary Table 8.xlsx"),sheet = i, skip = 3 )
  tmp <- tmp[,c(2,4)]
  colnames(tmp) <- c("Gene","Change")
  tmp <- tmp[tmp$Gene != "NA",]
  tmp <- as.data.frame(tmp[complete.cases(tmp),])
  row.names(tmp) <- tmp$Gene
  tmp$Source <- rep(paste0("Smajić.2022.", deg_names[k]),nrow(tmp))
  smajic_traj[[k]] <- tmp
 
}
names(smajic_traj) <- paste0("Smajić.2022.", deg_names)

# total
deg_list <- c(ferraro_deg_list,smajic_deg,smajic_traj)

# save marker genes
# save(deg_list, file=paste0(git_repo_location,"deg_list.published.R"))
