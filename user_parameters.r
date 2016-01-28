
# this script contains user-specific parameters referenced by the analysis scripts
# it is in .gitignore, so your local changes will not be pushed

# directories
project_stem_dir <- "C:/Lab_projects/2016_Phylomodelling"
climate_data_dir <- "E:/BCM/CA_2014/Summary/HST/Normals_30years"
filled_climate_data_dir <- "C:/Lab_projects/2016_Phylomodelling/Data/Climate/BCM_normals"
#occurrence_data_dir <- "E:/Phylo_modelling/Data/Species"
occurrence_data_dir <- paste0(project_stem_dir, "/Data/Species")
occurrence_data_dir_processed <- paste0(occurrence_data_dir, '/Processed3/')
maxent_output_dir <- "C:/Lab_projects/2016_Phylomodelling/Output/Maxent/"
#charts_output_dir <- "E:/Phylo_modelling/Output/Charts"

# number of nodes to use for parallel processing jobs
# (shoud be no greater than number of processors in your machine)
nodes <- 6

# Version of maxent run
maxent_run_version <- "V4"

# species list file name
#species_list <- paste0(spdir, '/0_Species_list_v2.rdata')

# occurrence background file name for maxent
#maxent_background <- paste(spdir, '10000_CA_plants_bg_810m.rdata', sep='')

# minimum number of records for maxent modeling
richness_min_cells <- 10

# maxent threshold statistic
maxent_threshold_stat <- "spec_sens"


