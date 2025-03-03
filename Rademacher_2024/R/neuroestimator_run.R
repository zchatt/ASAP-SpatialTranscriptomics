# NEUROeSTIMator - run Nakamura
# refer to neuroestimator_install.sh for install instructions

# Operating System: Ubuntu 22.04.4 LTS              
#           Kernel: Linux 6.8.0-48-generic
#     Architecture: x86-64

# # start R 
# cd /data/zac/neuroestimator
# conda activate neuroestimator
# R

# libs
library(remotes)
library(neuroestimator)

# #------------- run select spots from each mouse model ------------
# setwd("/data/zac/neuroestimator/select_spots")
# 
# count_files <- list.files(pattern=".rds")
# res <- list()
# 
# run_index <- 1:length(count_files)
# run_index <- run_index[!run_index %in% c(9,14,21)]
# 
# for (i in run_index) {
#   print(i)
#   counts <- readRDS(count_files[i])
#   res[[i]] <- neuroestimator(counts, species = "mmusculus")
#   #count_labels <- unlist(strsplit(gsub("transcripts", "", count_files[i]), "_")[[1]])
#   
# }
# 
# save(res,file="neuroestimator_res_291124.Rdata")

#------------- run all spots from viral DREADD model ------------
#setwd("/Users/zacc/USyd/spatial_transcriptomics/analysis/neuroestimator/250110_ForZac")
setwd("/data/zac/neuroestimator/250110_ForZac")

count_files <- list.files(pattern=".rds")
res <- list()

run_index <- 1:length(count_files)

for (i in run_index) {
  print(i)
  counts <- readRDS(count_files[i])
  #print(head(counts))
  #print(dim(counts))
  res[[i]] <- neuroestimator(counts, species = "mmusculus")
  #count_labels <- unlist(strsplit(gsub("transcripts", "", count_files[i]), "_")[[1]])
  
}

save(res,file="neuroestimator_res_130125.Rdata")
