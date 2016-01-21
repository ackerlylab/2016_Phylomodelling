
# this script contains user-specific parameters referenced by the analysis scripts
# it is in .gitignore, so your local changes will not be pushed

# directories
project_stem_dir <- "E:/Phylo_modelling/"
climate_data_dir <- "E:/BCM/CA_2014/Summary/HST/Normals_30years/"
occurrence_data_dir <- "E:/Phylo_modelling/Data/Species"

# number of nodes to use for parallel processing jobs
# (shoud be no greater than number of processors in your machine)
nodes <- 6

# minimum number of records for maxent modeling
maxent_min_records <- 10