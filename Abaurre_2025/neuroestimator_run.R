# NEUROeSTIMator - run Abaurre reference
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

setwd("/data/zac/neuroestimator")

#------------ single-cell ------------#
# read in count data
# counts <- read.table("/data/zac/neuroestimator/abaurre_countmmatrix_121224.txt")
# counts <- read.table("/data/zac/neuroestimator/20250122-abaurre_countmmatrix.txt")
counts <- read.table("/data/zac/neuroestimator/merged_countmatrix_DA_130224.txt")

# run neuroestimator
res <- neuroestimator(counts, species = "hsapiens")

# save
#save(res,file="neuroestimator_res_abaurre_121224.Rdata")
#save(res,file="neuroestimator_res_20250122-abaurre.Rdata")
save(res,file="neuroestimator_res_20250213-abaurre.Rdata")


#------------ Spatial ------------#
# read in count data
counts <- read.table("/data/zac/neuroestimator/V52Y16_079_A1_countmmatrix_121224.txt")

# run neuroestimator
res <- neuroestimator(counts, species = "hsapiens")

# save
save(res,file="neuroestimator_res_V52Y16_079_A1_121224.Rdata")
