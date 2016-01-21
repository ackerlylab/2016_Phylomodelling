

# Threshold Maxent predictions and sum to richness
# Matthew Kling
# January 2016


library(raster)
library(dismo)
library(doParallel)


######### file paths ##########
# source user_parameters.r before running #

clim_dir <- climate_data_dir
model_dir <- paste0(project_stem_dir, "/Output/Maxent/v3")
richness_dir <- paste0(project_stem_dir, "/Output/Richness")
pr_dir <- paste0(project_stem_dir, "/Data/Richness/point_richness_rasters")



######### climate data ###########
env.files <- list.files(path=clim_dir, pattern='.data', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt")) #use sort to ensure correct names used for predictors
files <- env.files[which(substr(env.files,5,8)%in%'810m')]


######### maxent models ###############
spp_dirs <- list.dirs(model_dir, full.names=T, recursive=F)


######### thresholded predictions for each species #############

predictBinary <- function(model, predictors){
        # function to determine threshold, make model prediction, and return threshold prediction
        eval <- evaluate(model@presence, model@absence, model)
        thresh <- threshold(eval, stat="spec_sens")
        pred <- predict(m, predictors)
        pred <- reclassify(pred, c(0, thresh, 0, thresh, 1, 1))
        return(pred)
}

cl <- makeCluster(nodes)
registerDoParallel(cl)

results <- foreach(spp = spp_dirs,
                   .packages=c("raster", "dismo")) %dopar% {
                           
                           m <- readRDS(paste0(spp, "/ModelObject.rdata"))
                           if(class(m)=="try-error") return("no maxent model")
                           p <- stack(lapply(files,function(x) readRDS(paste(clim_dir,x,sep=""))))
                           names(p) <- climnames
                           bp <- predictBinary(m, p)
                           saveRDS(bp, paste0(spp, "/BinaryRangePrediction.rds"))
                           return("success")
                   }

stopCluster(cl)


######## richness map summed across species ###########

# some spp failed due to nonexistent models; skip them
good <- file.exists(paste0(spp_dirs, "/BinaryRangePrediction.rds"))
spp_dirs <- spp_dirs[good]

sumRasters <- function(dir, pattern){
        # function returns the sum of all raster files in dir (and subdirs) that match pattern
        files <- list.files(dir, pattern=pattern, full.names=T, recursive=T)
        richness <- readRDS(files[1])
        for(f in files[2:length(files)]){
                f <- readRDS(f)
                if(class(f)!="RasterLayer") next()
                richness <- richness + f
        }
        return(richness)
}

richness <- sumRasters(model_dir, "BinaryRangePrediction.rds")
writeRaster(richness, paste0(richness_dir, "/richness.tif"), overwrite=T)


######### plot #########

library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)

r <- as.data.frame(rasterToPoints(richness))

# map
p <- ggplot(r, aes(x, y, fill=layer)) +
        geom_raster() +
        scale_fill_viridis() +
        coord_fixed(ratio=1) +
        theme(panel.background=element_blank(), panel.grid=element_blank(),
              axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
              legend.position="none", title=element_text(size=25)) +
        guides(fill = guide_colorbar(barwidth=20)) +
        labs(title=paste0("Species richness (", length(spp_dirs), " MaxEnt models)"))

# histogram data
h <- r %>%
        mutate(lyr = plyr::round_any(layer, 25)) %>%
        group_by(lyr) %>%
        summarize(n=n())

# histogram
l <- ggplot(h, aes(lyr, n, fill=lyr)) + 
        geom_bar(stat="identity", width=25) + 
        scale_fill_viridis() +
        theme(panel.background=element_blank(), panel.grid=element_blank(),
              axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none") +
        labs(x="number of species")

png(paste0(richness_dir, "/richness.png"), width = 1000, height = 1500)
plot(p)
print(l, 
      vp=viewport(x = .5, y = .55, 
                  width = unit(0.45, "npc"), height = unit(0.3, "npc"),
                  just = c("left", "bottom")))
dev.off()





############### comparison to occurrence density from David ###############


### create 25k and 50k upscaled versions of binary range maps ###

upscale <- function(file, template, tag){
        outfile <- sub("\\.rds", paste0("_", tag, ".rds"), file)
        if(file.exists(outfile)) return("file already exists")
        f <- readRDS(file)
        f <- as.data.frame(rasterToPoints(f))
        coordinates(f) <- c("x", "y")
        f <- rasterize(f, template, field="layer", fun=function(x, ...){max(na.omit(x))})
        saveRDS(f, outfile)
        return("upscaled")
}

pr_files <- list.files(pr_dir, pattern=".grd", full.names=T)

cl <- makeCluster(7)
registerDoParallel(cl)
results <- foreach(spp = paste0(spp_dirs, "/BinaryRangePrediction.rds"),
                   .packages=c("raster")) %dopar% {
                           if(!file.exists(spp)) return("no maxent model")
                           upscale(spp, template=raster(pr_files[1]), "25k")
                           upscale(spp, template=raster(pr_files[2]), "50k")
                           return("success")
                   }
stopCluster(cl)


### sum to create richness maps at both resolutions ###

writeRaster(sumRasters(model_dir, "BinaryRangePrediction_25k.rds"), 
            paste0(richness_dir, "/richness_25k.tif"), overwrite=T)
writeRaster(sumRasters(model_dir, "BinaryRangePrediction_50k.rds"),
            paste0(richness_dir, "/richness_50k.tif"), overwrite=T)


### compare to occurrence density maps ###

for(res in c("25", "50")){
        r <- stack(raster(paste0(richness_dir, "/richness_", res, "k.tif")),
                     stack(pr_files[grepl(paste0(res, "k"), pr_files)]))
        names(r) <- c("maxent", "rarefied", "raw")
        r <- mask(r, r[[3]])
        #r <- mask(r, calc(r, sum))
        r <- as.data.frame(rasterToPoints(r))
        
        s <- r %>%
                gather(stat, value, -x, -y) 
        
        p1 <- ggplot(s, aes(x, y, fill=value)) + 
                geom_raster() + 
                facet_grid(.~stat) +
                scale_fill_viridis(na.value="white")  +
                coord_fixed(ratio=1) +
                theme(panel.background=element_blank(), panel.grid=element_blank(),
                      axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                      legend.position="top", title=element_text(size=25)) +
                guides(fill = guide_colorbar(barwidth=20)) +
                labs(fill="species richness  ")
        
        t <- s %>%
                group_by(stat) %>%
                mutate(value = value/max(value, na.rm=T))
        
        p2 <- ggplot(t, aes(x, y, fill=value)) + 
                geom_raster() + 
                facet_grid(.~stat) +
                scale_fill_viridis(na.value="white")  +
                coord_fixed(ratio=1) +
                theme(panel.background=element_blank(), panel.grid=element_blank(),
                      axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                      legend.position="top", title=element_text(size=25)) +
                guides(fill = guide_colorbar(barwidth=20)) +
                labs(fill="normalized richness (% of max)  ")
        
        u <- t %>%
                spread(stat, value) %>%
                mutate(maxent=maxent-maxent, rarefied=rarefied-maxent, raw=raw-maxent) %>%
                gather(stat, value, -x, -y)
        
        p3 <- ggplot(u, aes(x, y, fill=value)) + 
                geom_raster() + 
                facet_grid(.~stat) +
                scale_fill_gradientn(colours=c("darkblue", "dodgerblue", "gray", "orange", "darkred"),
                                     limits=c(-1, 1), na.value="white")  +
                coord_fixed(ratio=1) +
                theme(panel.background=element_blank(), panel.grid=element_blank(),
                      axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                      legend.position="top", title=element_text(size=25)) +
                guides(fill = guide_colorbar(barwidth=20)) +
                labs(fill="difference from maxent (normalized)  ")
        
        p <- arrangeGrob(p1, p2, p3, ncol=1)
        
        png(paste0(paste0(richness_dir, "/richness_comparison_", res, ".png")), width=1000, height=1500)
        grid.draw(p)
        dev.off()
        
        png(paste0(paste0(richness_dir, "/richness_scatterplots_", res, ".png")), width=500, height=500)
        pairs(r[c("maxent", "rarefied", "raw")])
        dev.off()
}










######### species range sizes ############

cl <- makeCluster(7)
registerDoParallel(cl)
r <- foreach(spp = spp_dirs,
                   .packages=c("raster")) %dopar% {
                           if(!file.exists(paste0(spp, "/BinaryRangePrediction.rds"))) return("no data")
                           rasters <- paste0(spp, "/BinaryRangePrediction", c("", "_25k", "_50k"), ".rds")
                           ranges <- sapply(rasters, function(x) sum(na.omit(values(readRDS(x)))))
                           return(ranges)
                   }
stopCluster(cl)

good <- sapply(r, length)==3
rd <- r[good]
rd <- lapply(rd, as.vector)
rd <- do.call("rbind", rd)
rd <- cbind(basename(spp_dirs[good]), as.data.frame(rd))
names(rd) <- c("species", "range_810m", "range_25km", "range_50km")
write.csv(rd, paste0(richness_dir, "/range_size.csv"), row.names=F)
