# Validation of DEG within substantia nigra of parkinsons disease

library(readxl)

# repo location
git_repo_location = "/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/deg_markers/tables/"

#1) Ferraro et al., Cells, 2022
# A cause of the limited overlap in DEGs from published microarray studies (n=9) is the
# heterogeneity in the neurodegeneration and technical artefacts. 
# Meta-analysis of SNpc microarray datasets controlling for cell-type
# report LM results with a cell-type aware and unaware model
tmp <- read_excel(paste0(git_repo_location,"CorrectingForCellPropPDsn_SuppTableS3.xlsx"),sheet = 1, skip = 2 )

tmp2 <- tmp[tmp$is_cell_unaware_LMM_DEG == "yes",]
tmp2 <- tmp2[,c("Gene","cell_unaware_LMM_estimate")]
colnames(tmp2) <- c("Gene","Change")
tmp2$Source <- rep("Ferraro.2022.cellunaware",nrow(tmp2))
ferraro_deg_cellunaware <- tmp2

tmp2 <- tmp[tmp$is_cell_aware_LMM_DEG == "yes" ,]
tmp2 <- tmp2[,c("Gene","cell_aware_LMM_estimate")]
colnames(tmp2) <- c("Gene","Change")
tmp2$Source <- rep("Ferraro.2022.cellaware",nrow(tmp2))
ferraro_deg_cellaware <- tmp2

#2) Smajić et al., Brain, 2022
# snRNAseq of PD and CTR SNpc
# ID population of DaN with high CADPS2
# increased number of activated microglia
# increased number of astrogliosis

# Differential up-regulated genes for cell-types trajectories and in PD for each cell-type
deg_names <- c(#"Microglia_trajectory_genes",
               "Microglia_IPD_diff_exp_genes",
               #"Astrocyte_trajectory_genes",
               "Astrocyte_IPD_diff_exp_genes",
               #"Oligodendrocyte_trajectory_genes",
                "Oligodendrocytes_IPD_diff_exp_genes",
                "Inhibitory_IPD_diff_exp_genes",
                "Excitatory_IPD_diff_exp_genes",
                "GABA_IPD_diff_exp_genes",
                "DaNs_IPD_diff_exp_genes",
                "CADPS2high_IPD_diff_exp_genes",
                "OPC_IPD_diff_exp_genes",
                "Ependymal_IPD_diff_exp_genes",
                "Pericytes_IPD_diff_exp_genes",
                "Endothelial_cell_IPD_diff_exp_genes")
smajic_deg <- list()
for (i in c(2,4,6:15)){
  tmp <- read_excel(paste0(git_repo_location,"Supplementary Table 8.xlsx"),sheet = i, skip = 3 )
  tmp <- tmp[,c(1,4)]
  colnames(tmp) <- c("Gene","Change")
  tmp$Source <- rep(paste0("Smajić.2022.", deg_names[i]),nrow(tmp))
  smajic_deg[[i]] <- tmp
}

deg_names <- c("Microglia_trajectory_genes",
               "Microglia_IPD_diff_exp_genes",
               "Astrocyte_trajectory_genes",
               "Astrocyte_IPD_diff_exp_genes",
               "Oligodendrocyte_trajectory_genes",
               "Oligodendrocytes_IPD_diff_exp_genes",
               "Inhibitory_IPD_diff_exp_genes",
               "Excitatory_IPD_diff_exp_genes",
               "GABA_IPD_diff_exp_genes",
               "DaNs_IPD_diff_exp_genes",
               "CADPS2high_IPD_diff_exp_genes",
               "OPC_IPD_diff_exp_genes",
               "Ependymal_IPD_diff_exp_genes",
               "Pericytes_IPD_diff_exp_genes",
               "Endothelial_cell_IPD_diff_exp_genes")
smajic_deg <- list()
for (i in 1:15){
  tmp <- read_excel(paste0(git_repo_location,"Supplementary Table 8.xlsx"),sheet = i, skip = 3 )
  tmp <- tmp[,c("Gene name","estimate")]
  colnames(tmp) <- c("Gene","Change")
  tmp$Source <- rep(paste0("Smajić.2022.", deg_names[i]),nrow(tmp))
  smajic_deg[[i]] <- tmp
}


# total
deg_df <- rbind(ferraro_deg_cellunaware,ferraro_deg_cellaware,smajic_deg)

# save marker genes
# save(deg_list_total, file=paste0(git_repo_location,"deg_list.published.R"))
