# Quality assessment of Visium human probes in mouse tissue for running of Xenografts on Visium 10X human array.

# Libraries
library(stringr)
library(DECIPHER)
library(seqinr)
library(Rsamtools)
library(ggplot2)
libraru(biomaRt)

#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}


# setwd
setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/xenograft")

#####################
### Pre-alignment ###
#####################

# read in probe sets
human_v1 <- read.delim("/Users/zacc/USyd/spatial_transcriptomics/data/probe_sets/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv",
                       sep=",", skip=5)
human_v1$gene <- str_match(human_v1$probe_id, "\\|(.*?)\\|")[,2]

human_v2 <- read.delim("/Users/zacc/USyd/spatial_transcriptomics/data/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv",
                       sep=",", skip=5)
human_v2$gene <- str_match(human_v2$probe_id, "\\|(.*?)\\|")[,2]

mouse_v1 <- read.delim("/Users/zacc/USyd/spatial_transcriptomics/data/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv",
                       sep=",", skip=5)
mouse_v1$gene <- str_match(mouse_v1$probe_id, "\\|(.*?)\\|")[,2]

# # convert mouse gene names to human equivalent. Taken from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/ - Thank you!
## Ensemble servers down  - will revisit
# musGenes <- mouse_v1$gene
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="asia.ensembl.org")
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="asia.ensembl.org")
# genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = musGenes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
# humanx <- unique(genesV2[, 2])

# we are using the human v2 probe set
probe_set <- read.delim("/Users/zacc/USyd/spatial_transcriptomics/data/SANPIN_VisiumFFPE_Cytassist_results/221011/V52Y16-079-A1/VISIUM/V52Y16-079-A1/outs/probe_set.csv",
                       sep=",", skip=5)

# probe overlap between human and mouse
stat_collect <- matrix(0,10,2)
stat_collect[1,] <- table(mouse_v1$probe_seq %in% human_v2$probe_seq)
stat_collect[2,] <- table(human_v2$probe_seq %in% mouse_v1$probe_seq)

# write probes to FASTA file for alignment with STAR
write.fasta(as.list(human_v2$probe_seq), human_v2$probe_id, file.out = "human_v2_probes.fasta", open = "w", nbchar = 60, as.string = FALSE)
write.fasta(as.list(mouse_v1$probe_seq), mouse_v1$probe_id, file.out = "mouse_v1_probes.fasta", open = "w", nbchar = 60, as.string = FALSE)

# Please refer to "align_species.pbs" script for details on alignment 

#####################
### Post-alignment ###
#####################
# read in aligned BAM files
bamFile <- BamFile("human_map_humanAligned.sortedByCoord.out.bam")
hh_bam <- scanBam(bamFile)
bamFile2 <- BamFile("mouse_map_humanAligned.sortedByCoord.out.bam")
mh_bam  <- scanBam(bamFile2)

# MAPQ distribution
MAPQ <- as.numeric(c(hh_bam[[1]]$mapq, mh_bam[[1]]$mapq))
mapped_probes <- c(rep("human:human",length(hh_bam[[1]]$mapq)),rep("mouse:human",length(mh_bam[[1]]$mapq)))

df <- as.data.frame(cbind(mapped_probes,MAPQ))
plot_hist_mapq <- ggplot(df, aes(x=MAPQ, color=mapped_probes)) +
  geom_histogram(fill="white", alpha=0.5, stat="count",position="dodge") + theme_minimal() +
  stat_count(binwidth = 1, 
             geom = 'text', 
             aes(label = ..count..),
             position = position_stack(vjust = 0.8))

# Probes muniquely apping mouse genome
#store names of BAM fields
bam_field <- names(mh_bam[[1]])
#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(mh_bam, "[[", y)))
#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field
# select unique mapping
mh_bam_df <- bam_df
mh_bam_df <- mh_bam_df[!is.na(mh_bam_df$mapq),]
mh_bam_df <- mh_bam_df[mh_bam_df$mapq == 255,]
# get gene names
mh_bam_df$gene <- str_match(mh_bam_df$qname, "\\|(.*?)\\|")[,2]

# number of unique genes mapped = 4731
length(unique(mh_bam_df$gene))

# number of genes with same name as human = 4340 T, 391 F
table(unique(toupper(mh_bam_df$gene)) %in% human_v2$gene)

# number of probes mapping same genes as human = 1 T, 4778 F
h_gs <- paste0(human_v2$probe_seq,sep=".",human_v2$gene)
m_gs <- paste0(mh_bam_df$seq,sep=".",toupper(mh_bam_df$gene))
table(m_gs %in% h_gs)

# load marker gene lists
git_repo_location = "/Users/zacc/github_repo/ASAP-SpatialTranscriptomics/celltype_markers/tables/"
load(paste0(git_repo_location,"gene_marker_list.cellenrich.R"))

# calculate numbers of cell-type marker genes covered by human probes that uniquely map mouse genome



# read in BED files of reference FASTA matching aligned BAM
hh_bed <- read.delim(file="human_map_humanAligned.sortedByCoord.out.bed", sep="\t", header=FALSE)

# quantify delta-G for each alignment of human 


# delta-G cutoff for each unique alignment v else human


# quantify delta-G for each unique alignment of mouse


# hamming distance matrix 
myXStringSet<- DNAStringSet(c(mouse_v1$probe_seq,human_v2$probe_seq))
dist_mat <- DistanceMatrix(myXStringSet,correction = "none", method = "longest")
# multiply by probe length to get hamming distance
dist_mat <- dist_mat * 



# NOTES
# we could identify a good list of probes mapping to the mouse genome through alignment to 
# rat genome (bwa-mem) and then delta-G calculation (primer-3)




# list of probes to keep



# Filter the human v2 probes to create a mouse applicable CSV
# Eg. Edit the probes CSV file so the probe’s 'included' column is 'TRUE'. This selectively includes the probe in spaceranger count analysis.



# Xenograft samples will be aligned twice, once using standard human V2 probeset and again using mouse applicable probe set

spaceranger count \
  --id=V52Y16-079-B1 \
  --transcriptome=/directflow/GWCCGPipeline/projects/reference/refdata-cellranger-GRCh38-2020-A \
  --fastqs=/directflow/GWCCGPipeline/projects/bioinformatics/SANPIN_VisiumFFPE_Cytassist/fastq_path/221007_A00152_0716_AH7KLFDSX5/V52Y16-079-B1/VISIUM \
  --sample=V52Y16-079-B1 \
  --slide=V52Y16-079 \
  --area=B1 \
  --image=resources/count/SANPIN_VisiumFFPE_Cytassist/V52Y16-079-B1/VISIUM/V52Y16-079-B1.tif \
  --cytaimage=resources/count/SANPIN_VisiumFFPE_Cytassist/V52Y16-079-B1/VISIUM/V52Y16-079-B1.ca.tif \
  --loupe-alignment=resources/count/SANPIN_VisiumFFPE_Cytassist/V52Y16-079-B1/VISIUM/V52Y16-079-B1.json \
  --probe-set=resources/count/SANPIN_VisiumFFPE_Cytassist/V52Y16-079-B1/VISIUM/probe_set.csv \
  --jobmode=local \
  --localcores=16 \
  --localmem=256

spaceranger count \
--id=V52Y16-079-B1 \
--transcriptome=/directflow/GWCCGPipeline/projects/reference/refdata-cellranger-GRCh38-2020-A \
--fastqs=/directflow/GWCCGPipeline/projects/bioinformatics/SANPIN_VisiumFFPE_Cytassist/fastq_path/221007_A00152_0716_AH7KLFDSX5/V52Y16-079-B1/VISIUM \
--sample=V52Y16-079-B1 \
--slide=V52Y16-079 \
--area=B1 \
--image=resources/count/SANPIN_VisiumFFPE_Cytassist/V52Y16-079-B1/VISIUM/V52Y16-079-B1.tif \
--cytaimage=resources/count/SANPIN_VisiumFFPE_Cytassist/V52Y16-079-B1/VISIUM/V52Y16-079-B1.ca.tif \
--loupe-alignment=resources/count/SANPIN_VisiumFFPE_Cytassist/V52Y16-079-B1/VISIUM/V52Y16-079-B1.json \
--probe-set=Visium_Human_mouseapplicable_Transcriptome_Probe_Set_v2.0_GRCh38-2023-A.csv \
--jobmode=local \
--localcores=16 \
--localmem=256