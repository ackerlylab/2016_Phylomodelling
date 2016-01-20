# Script for running simple MaxEnt models for California plant
# Output includes maxent modelling files, and .rdata file for each species
# including list with training AUC, maximum suitability value, suitability
# values for each occurrence ID and ratio suitability/max suitability


# Clear workspace
rm(list=ls())
start=Sys.time()

#-----------------#
# Set directories #
#-----------------#
Computer <- 'Ackerly'

if(Computer == 'Ackerly') {
        wdir <- 'E:/Phylo_modelling/'
        cdir <- 'E:/BCM/CA_2014/Summary/HST/Normals_30years/'
} else if (Computer == 'HP') {
        wdir <- 'D:/Phylo_modelling/'
        cdir <- 'D:/BCM/CA_2014/Summary/HST/Normals_30years/'
}

spdir <- paste(wdir, 'Data/Species/Processed/', sep='')
#bgfile <- paste(spdir, '10000_CA_plants_bg.rdata', sep='')
bgfile <- paste(spdir, '10000_CA_plants_bg_1080m.rdata', sep='')
odir <- paste(wdir, 'Output/', sep='')

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

# Arguments for maxent models
mxArgs <- c("-a", "-z", "outputformat=raw", "maximumbackground=10000", 
            "nothreshold", "nohinge")
# What species?
allSpecies <- readRDS(paste(spdir, '0_Species_list.rdata', sep=''))

#--------------------------------------------#
# Run full models with filled climate layers #
#--------------------------------------------#
# Updated: running on 1080 m resolution instead


#files <- env.files[which(substr(env.files,5,10)%in%'filled')]
files <- env.files[which(substr(env.files,5,9)%in%'1080m')]
predictors <- stack(lapply(files,function(x) readRDS(paste(cdir,x,sep=""))))
names(predictors) = climnames 
orig.project <- '+proj=longlat +ellps=WGS84'
bg <- readRDS(bgfile)

library(doParallel)
cl <- makeCluster(7)
registerDoParallel(cl)

results <- foreach(i=1:length(allSpecies)) %dopar% {
#for (i in 1:length(maxent_failed)) {
        options(java.parameters = "-Xmx1g" )
        require(dismo)
        require(rJava)
        require(raster) 
        extract <- raster::extract
        
        mySpecies <- allSpecies[i]
        #mySpecies <- maxent_failed[i]
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
        #mx.dir = paste(odir, 'Maxent/V1/',mySpecies, sep="")
        mx.dir = paste(odir, 'Maxent/V2/',mySpecies, sep="")
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

