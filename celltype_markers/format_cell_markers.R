### Published gene marker sets that we use for cell-type enrichment (PAGE) of spatial transcriptomic data
library(readxl)
library(biomaRt)

# repo location
git_repo_location = "/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/celltype_markers/tables/"

## 0) cell-type marker database
# # NOTE: Adrenergic markers were limited. Also database doesn't say where markers were originally identified.
# library(clustermole)
# markers <- clustermole_markers(species = "hs")
# tmp <- markers[markers$organ == "Brain",]
# tmp$gene_original[tmp$celltype == "Adrenergic neurons"] %in% webber_gl$gene_name
# # NOTE: Adrenergic markers were limited. Also database doesn't say where markers were originally identified.

## 1) Gene lists of differential expression analysis for ten dopaminergic subtypes (snRNAseq) reported in Kamath et al, 2022 - https://pubmed.ncbi.nlm.nih.gov/35513515/
kamath_gl <- read_excel(paste0(git_repo_location,"41593_2022_1061_MOESM3_ESM.xlsx"),sheet = "Supplementary_Table_8" )
kamath_gl <- as.data.frame(kamath_gl)
kamath_gl <- kamath_gl[,c("DA_subtype","primerid","fdr")]
colnames(kamath_gl) <- c("cell","gene_name","FDR")
kamath_gl$study <- rep("Kamath.et.al.2022",nrow(kamath_gl))

## 2) Gene lists from single cell analysis of Norepinephrine (NE) neurons in the locus coeruleus (LC) from Weber et al, 2022 - https://elifesciences.org/reviewed-preprints/84628
# 327 statistically significant genes with elevated expression in the NE neuron cluster. 
webber_gl <- read.table(paste0(git_repo_location,"media-4.csv"),sep=",", header=TRUE)
webber_gl <- as.data.frame(webber_gl[webber_gl$FDR < 0.05,])
webber_gl <- cbind(cell=rep("NE",nrow(webber_gl)),webber_gl[,c("gene_name","FDR")])
webber_gl <- webber_gl[!webber_gl$gene_name %in% kamath_gl$gene_name,]
webber_gl$study <- rep("Webber.et.al.2022",nrow(webber_gl))

## 3) Gene lists 10 substantia nigra cell populations identified by snRNAseq from Agarwal et al, 2020 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7442652/
agarwal_gl <- read_excel(paste0(git_repo_location,"41467_2020_17876_MOESM6_ESM.xlsx"),sheet = "format_zc" )
agarwal_gl <- as.data.frame(agarwal_gl )
agarwal_gl <- agarwal_gl[,c("Cluster","Gene","False Discovery Rate" )]
colnames(agarwal_gl) <- c("cell","gene_name","FDR")
agarwal_gl$study <- rep("Agarwal.et.al.2020",nrow(agarwal_gl))

## 4) Gene lists from meta-analysis of tumor-infiltrating T-cells defining conserved human T-cell subtypes (n=9) gene sets from Andreatta et al, 2021 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8137700/
andreatta_gl <- read_excel(paste0(git_repo_location,"41467_2021_23324_MOESM4_ESM.xlsx"),sheet = "format_zc" )
andreatta_gl <- as.data.frame(andreatta_gl)
colnames(andreatta_gl) <- c("gene_name","cell")
andreatta_gl$FDR <- rep(NA,nrow(andreatta_gl))
andreatta_gl$study <- rep("Andreatta.et.al.2021",nrow(andreatta_gl))

## 5) Gene lists from snRNA-seq of white matter of 6 oligodendrocyte subtypes from Jäkel et al, 2019 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6544546/
olig_names <- c("OPC","COP","Oligo1","Oligo2","Oligo3","Oligo4","Oligo5","Oligo6","imOLG")
ol_list <- list()
for (i in 2:10){
  tmp <- read_excel(paste0(git_repo_location,"EMS81327-supplement-Supplementary_Table_4.xlsx"),sheet = i )
  tmp$cell_type <- rep(olig_names[i-1],nrow(tmp))
  ol_list[[i-1]] <- tmp
}

jakel_gl <- as.data.frame(do.call("rbind", ol_list))
jakel_gl <- jakel_gl[,c(6,7)]
colnames(jakel_gl ) <- c("gene_name","cell")
jakel_gl$FDR <- rep(NA,nrow(jakel_gl))
jakel_gl$study <- rep("Jakel.et.al.2019",nrow(jakel_gl))

# 6) Gene lists manually curated for iPSC derived neuron scRNAseq cell-type assignment used in Jerber et al, 2021 - https://www.nature.com/articles/s41588-021-00801-6
Proliferating = c("MKI67", "TOP2A", "KIAA1524")
Neural_Prog = c("LIX1", "RAX", "NR2F1", "NES")
FP_prog = c("ZEB2", "DMBX1","HMGA1", "HMGB2")
Glia = c("TNC", "SOX2", "CDH2", "HES1")
FP = c("LMX1A", "FOXA2")
Astrocytes = c("S100B", "AQP4", "GFAP", "SLC1A3","SOX9")
Ependyma = c("STOML3","CCDC153","CDHR4","FOXJ1","DNAH11","TTR","MLF1")
Serotonergic = c("CHGB", "DDC","FEV","GATA2", "GATA3", "GCH1", "GCHFR", "HTR1A","HTR1B","MAOA","MAOB","SLC18A2","SLC29A4","SLC6A4","TPH2")
Neuron = c("SYT1", "SNAP25")
Neuroblasts = c("NEUROG1", "NEUROD1", "NEUROG2", "NHLH1", "SIM1")
Dopaminergic =c("ABCC8","ACOT7","ALDH1A1","AMER3","ARG2","ASB4","BNC2","CADPS2","CALB1","CALB2","CAMK2N1",
                            "CCK","CDK14","CDKN1C","CHL1",
                            "CHRNA4","CPEB3","CXCR4","DCC","DKK3","DRD2",
                            "EBF2","EN1","EN2","EPHA5","ERC2","FGF13","FOXA1","FOXA2","GDAP1","GFRA1","GRIA3","GRIK3",
                            "GRP","ICA1L","IGF1","KCNIP4","KCNJ6","KLHL1","KLHL13","LGI1","LMO3","LMX1A","LMX1B",
                            "LRRC3B","LRRTM2","LSAMP","LY6H","NETO2","NR4A2","NTSR1","OTX2","PBX1","PRICKLE2","PRKCA",
                            "PRL","PRRT4","PITX3","PTPN5","PTPRO","RET","SCG2","SLC10A4","SLC18A1","SLC18A2","SLC18A3",
                            "SLC6A3","SNCA","SOX6","TH","TMCC3","TMEFF2","TMEM255A","TUB","VGF")
Cortical_hem =c("DKK3","EOMES","RELN","LHX1","LHX5",
                            "TP73","CALB2","EBF3","SAMD3","SOX2")

gene_name <- c(Proliferating,Neural_Prog,FP_prog,Glia,FP,Astrocytes,Ependyma,Serotonergic,Neuron,Neuroblasts,Dopaminergic,Cortical_hem)
cell <- c(rep("Proliferating",length(Proliferating)),
  rep("Neural_Prog",length(Neural_Prog)),
  rep("FP_prog",length(FP_prog)),
  rep("Glia",length(Glia)),
  rep("FP",length(FP)),
  rep("Astrocytes",length(Astrocytes)),
  rep("Ependyma",length(Ependyma)),
  rep("Serotonergic",length(Serotonergic)),
  rep("Neuron",length(Neuron)),
  rep("Neuroblasts",length(Neuroblasts)),
  rep("Dopaminergic",length(Dopaminergic)),
  rep("PCortical_hem",length(Cortical_hem)))
FDR <- rep(NA,length(cell))
study <- rep("Jerber.et.al.2023",length(cell))

jerber_gl <- as.data.frame(cbind(gene_name,cell,FDR,study))

# 7) Top genes expressed in genes of 3 broad DaN subtypes in mouse from snRNAseq of SNpcGene reported in Azcorra et al, Nat Neuro, 2023 - https://www.nature.com/articles/s41593-023-01401-9
# Slc17a6 (Vglut2), Calb1 and Anxa1 cell-types
# Calb1+ and VGlut2+ subtypes; positive responses to unexpected rewards and aversive stimuli, deceleration correlated
# Distinct from Anxa1+ subtype;  acceleration correlated
# mouse IDs
Vglut2_DaN = c("Slc17a6", "Kcnip4", "Il1rapl2","Shisa6","Asic2")
Calb1_DaN = c("Calb1", "Fbxl7", "Cadm1","Pcsk2","Sdk1")
Anxa1_DaN = c("Anxa1", "Hs6st3", "Lmo3","Cntn5","Zfp804b")

# convert mouse to human IDs - Note; using archived ensembl due to getLDS issues - https://support.bioconductor.org/p/9144001/
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

musGenes <- Vglut2_DaN
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = musGenes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
Vglut2_DaN <- genesV2$HGNC.symbol

musGenes <- Calb1_DaN 
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = musGenes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
Calb1_DaN  <- genesV2$HGNC.symbol

musGenes <- Anxa1_DaN
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = musGenes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
Anxa1_DaN  <- genesV2$HGNC.symbol

gene_name <- c(Vglut2_DaN,Calb1_DaN,Anxa1_DaN)
cell <- c(rep("Vglut2_DaN",length(Vglut2_DaN)),
          rep("Calb1_DaN",length(Calb1_DaN)),
          rep("Anxa1_DaN",length(Anxa1_DaN)))
FDR <- rep(NA,length(cell))
study <- rep("Azcorra.et.al.2023",length(cell))
azcorra_gl <- as.data.frame(cbind(gene_name,cell,FDR,study))

###  combine into marker gene list and save
all_deg <- rbind(webber_gl,kamath_gl,agarwal_gl,andreatta_gl,jakel_gl,jerber_gl,azcorra_gl)
all_deg$cell_study <- paste0(all_deg$cell,sep=".",all_deg$study)
table(all_deg$cell_study)

# total
cells2 <- unique(all_deg$cell_study)
gene_marker_list <- list()
for (i in 1:length(cells2)){
  tmp <- all_deg[all_deg$cell_study %in% cells2[i],]
  gene_marker_list[[i]] <- tmp$gene_name
}
names(gene_marker_list) <- c(cells2)

# save marker genes
save(gene_marker_list, file=paste0(git_repo_location,"gene_marker_list.cellenrich.R"))
