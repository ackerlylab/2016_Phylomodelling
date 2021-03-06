
# this script contains user-specific parameters referenced by the analysis scripts
# it is in .gitignore, so your local changes will not be pushed

# directories

# For lab computer
project_stem_dir <- "E:/Phylo_modelling/"
climate_data_dir <- "E:/BCM/CA_2014/Summary/HST/Normals_30years/"
occurrence_data_dir <- paste0(project_stem_dir, "Data/Species")
occurrence_data_dir_processed <- paste0(occurrence_data_dir, '/Processed3/')
maxent_output_dir <- "E:/Phylo_modelling/Output/Maxent"

charts_output_dir <- "E:/Phylo_modelling/Output/Charts"

# Version of maxent run
maxent_run_version <- "V3"

# species list file name
species_list <- paste0(spdir, '0_Species_list_v2.rdata')

# occurrence background file name for maxent
maxent_background <- paste(spdir, '10000_CA_plants_bg_810m.rdata', sep='')

# number of nodes to use for parallel processing jobs
# (shoud be no greater than number of processors in your machine)
nodes <- 6

# minimum number of records for maxent modeling
maxent_min_records <- 10

