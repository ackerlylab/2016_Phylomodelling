
# this script contains user-specific parameters referenced by the analysis scripts
# it is in .gitignore, so your local changes will not be pushed

# directories
project_stem_dir <- "C:/Lab_projects/2016_Phylomodelling"
climate_data_dir <- "E:/BCM/CA_2014/Summary/HST/Normals_30years"
filled_climate_data_dir <- "C:/Lab_projects/Data/Climate/BCM_normals"
occurrence_data_dir <- "E:/Phylo_modelling/Data/Species"
maxent_output_dir <- "E:/Phylo_modelling/Output/Maxent"
charts_output_dir <- "E:/Phylo_modelling/Output/Charts"

# number of nodes to use for parallel processing jobs
# (shoud be no greater than number of processors in your machine)
nodes <- 6

# minimum number of records for maxent modeling
maxent_min_records <- 10