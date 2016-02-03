# Script for running simple MaxEnt models for California plant
# Output includes maxent modelling files, and .rdata file for each species

# By Naia Morueta-Holme
# Last updated by Matthew Kling Jan 27 2015

# Clear workspace
#rm(list=ls())
start=Sys.time()

#-----------------#
# Set directories #
#-----------------#

# source user_parameters.r before running

cdir <- filled_climate_data_dir
spdir <- occurrence_data_dir_processed
odir <- maxent_output_dir
species_list <- paste0(spdir, 'combined/0_Species_list_v2.rdata')
#maxent_background <- paste(spdir, '10000_CA_plants_bg_810m.rdata', sep='')
maxent_background <- paste(spdir, 'combined/10000_CA_plants_bg_810m_occbias.rdata', sep='')


#----------------#
# Load libraries #
#----------------#
options(java.parameters = "-Xmx1g" )
require(dismo)
require(rJava)
require(raster) 

#----------------#
# Set parameters #
#----------------#
# Set the environmental file names
files <- list.files(path=cdir, pattern='810m_filled3x', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt")) #use sort to ensure correct names used for predictors
mxModelType <- paste(climnames,collapse="-")
#searchString <- '810m'
#files <- env.files[which(substr(env.files,5,8)%in%searchString)]

# Arguments for maxent models
mxArgs <- c("-a", "-z", "outputformat=raw", "maximumbackground=10000", 
            "nothreshold", "nohinge")

# Version of maxent run (to assign output directory)
version <- maxent_run_version

# What species?
splist <- species_list
allSpecies <- readRDS(splist)

# Background file
bgfile <- maxent_background

# Number of nodes for parallel computing
no.nodes <- nodes

#-----------------#
# Run full models #
#-----------------#

predictors <- stack(lapply(files,function(x) readRDS(paste(cdir,x,sep="/"))))
names(predictors) = climnames 
orig.project <- '+proj=longlat +ellps=WGS84'
bg <- readRDS(bgfile)

library(doParallel)
cl <- makeCluster(no.nodes)
registerDoParallel(cl)

results <- foreach(i=1:length(allSpecies)) %dopar% {
        #for (i in 1:length(allSpecies)) {
        options(java.parameters = "-Xmx1g" )
        require(dismo)
        require(rJava)
        require(raster) 
        
        mySpecies <- allSpecies[i]
        print(mySpecies)
        Sys.time()
        
        # Read in species occurrence data and background
        pres <- readRDS(paste0(spdir, "/atomic/", mySpecies, ".rdata"))
        coordinates(pres) <- ~longitude + latitude
        projection(pres) <- CRS(orig.project)
        occur <- spTransform(pres, projection(predictors))
        
        # Directory to write files to
        mx.dir <- paste(odir, version, mySpecies, sep="/")
        if(file.exists(mx.dir) == F) {dir.create(mx.dir, recursive=T)}
        
        # Drop occurrences that fall in the water
        valid <- !is.na(extract(predictors[[1]], occur))
        occur <- occur[valid,]
        
        # Run the model!
        mx <- try(maxent(predictors, p=occur, a=bg ,progress='text', path=mx.dir, args=mxArgs))
        saveRDS(mx, file = paste(mx.dir, 'ModelObject.rdata', sep='/'))
}

stopCluster(cl)
Sys.time()-start

