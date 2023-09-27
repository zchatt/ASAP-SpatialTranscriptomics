# Instructions for specifying python and giotto environments
# Run on Macbook Pro (Apple M2 Max) with MacOS Ventura Version 13.5.

# set python paths
my_python_path = '/Users/zacc/opt/anaconda3/envs/giotto_env/bin/python'
Sys.setenv(RETICULATE_PYTHON = my_python_path)
RETICULATE_PYTHON=my_python_path

# load
library(devtools)

# Ensure Giotto Suite is installed
if(!"Giotto" %in% installed.packages()) {
  devtools::install_github("drieslab/Giotto@suite")
}
library(Giotto)

if(!"GiottoUtils" %in% installed.packages()) {
  devtools::install_github("drieslab/GiottoUtils")
}
library(GiottoUtils)

if(!"GiottoClass" %in% installed.packages()) {
  devtools::install_github("drieslab/GiottoClass")
}
library(GiottoClass)

# Ensure the Python environment for Giotto has been installed
genv_exists = checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to install the Giotto environment
  installGiottoEnvironment()
}

my_instructions = createGiottoInstructions(python_path = my_python_path)
 
# # load
# library(reticulate)
# library(Giotto)
# 
# use_python(my_python_path, required = T)
# 
# use_condaenv("py311") # This is a python 3.11.4 environment
# 
# # giotto environment
# installGiottoEnvironment(
#   packages_to_install = c("pandas==1.5.1", "networkx==2.8.8", "python-igraph==0.10.2",
#                           "leidenalg==0.9.0", "python-louvain==0.16", "python.app==1.4", "scikit-learn==1.1.3"),
#   #python_version = "3.11.4",
#   force_miniconda = FALSE,
#   force_environment = FALSE,
#   verbose = TRUE
# )