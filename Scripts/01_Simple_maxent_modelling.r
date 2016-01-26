# Script for running simple MaxEnt models for California plant
# Output includes maxent modelling files, and .rdata file for each species
# including list with training AUC, maximum suitability value, suitability
# values for each occurrence ID and ratio suitability/max suitability

# Last updated Jan 21 2015 by Naia Morueta-Holme

# Clear workspace
rm(list=ls())
start=Sys.time()

#-----------------#
# Set directories #
#-----------------#

# source user_parameters.r before running
# or user_parameters_maxent_cleaning.r to load original cleaning parameters

cdir <- climate_data_dir
spdir <- occurrence_data_dir_processed
odir <- maxent_output_dir

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
# Set the environmental file names (right now set to California BCM layers)
env.files <- list.files(path=cdir, pattern='.data', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt")) #use sort to ensure correct names used for predictors
mxModelType <- paste(climnames,collapse="-")
files <- env.files[which(substr(env.files,5,5+nchar(searchString-1))%in%searchString)]

# Arguments for maxent models
mxArgs <- c("-a", "-z", "outputformat=raw", "maximumbackground=10000", 
            "nothreshold", "nohinge")

# Version of maxent run (to assign output directory)
version <- maxent_run_version

# What species?
splist <- species_list
allSpecies <- readRDS(species_list)

# Background file
bgfile <- maxent_background

# Number of nodes for parallel computing
no.nodes <- nodes

#------------------#
# Run full models  #
#------------------#

predictors <- stack(lapply(files,function(x) readRDS(paste(cdir,x,sep=""))))
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
        extract <- raster::extract
        
        mySpecies <- allSpecies[i]
        writeLines(mySpecies)
        Sys.time()
        
        # Read in species occurrence data and background
        pres <- readRDS(paste(spdir, mySpecies, ".rdata", sep=""))
        coordinates(pres) <- ~longitude + latitude
        projection(pres) <- CRS(orig.project)
        occur <- spTransform(pres, projection(predictors))
        
        occur <- occur[!is.na(extract(predictors[[1]], occur)),]
        if(nrow(occur)==0){
                writeLines("           no occurrences overlapping climate data")
                next()}
        
        
        # Directory to write files to
        mx.dir = paste(odir, version, mySpecies, sep="/")
        if(file.exists(mx.dir) == F) {dir.create(mx.dir, recursive=T)} 
        
        # Run the model!
        mx <- try(maxent(predictors, p=occur, a=bg ,progress='text', path=mx.dir, args=mxArgs))
        saveRDS(mx, file = paste(mx.dir, 'ModelObject.rdata', sep='/'))
        
        
        # Predict suitability at occurrences
        cp <- extract(predictors, occur)
        
        predp <- predict(mx, cp)
        preda <- predict(mx, mx@absence)
        
        # Save results to .rdata object
        results <- list()
        results[['summary']] <- data.frame(species=mySpecies, auc.train=as.numeric(mx@results['Training.AUC',]),records.used=nrow(mx@presence))
        results[['occur']] <- occur
        results[['occur']]$suitability <- predp
        results[['occur']] <- cbind(results$occur, as.data.frame(cp))
        results[['summary']]$max.suit.occ <- max(predp, na.rm=T)
        results[['bg']] <- bg
        results[['bg']]$bgsuitability <- preda
        results[['summary']]$max.suit.all <- max(c(predp, preda), na.rm=T)
        results[['occur']]$prop_max.suit.occ <- results[['occur']]$suitability / results[['summary']]$max.suit.occ
        results[['occur']]$prop_max.suit.all <- results[['occur']]$suitability / results[['summary']]$max.suit.all
        
        saveRDS(results, paste(mx.dir, 'Suitability_data.rdata',sep='/'))
}

stopCluster(cl)
Sys.time()-start

