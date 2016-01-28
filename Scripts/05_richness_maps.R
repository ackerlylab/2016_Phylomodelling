

# Threshold Maxent predictions and sum to richness
# Matthew Kling
# January 2016


library(raster)
library(dismo)
library(doParallel)


######### file paths ##########
# source user_parameters.r before running #

clim_dir <- filled_climate_data_dir
model_dir <- paste0(project_stem_dir, "/Output/Maxent/v4")
richness_dir <- paste0(project_stem_dir, "/Output/Richness")
pr_dir <- paste0(project_stem_dir, "/Data/Richness/point_richness_rasters")


######### climate data ###########
files <- list.files(path=clim_dir, pattern='810m_filled3x', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt")) #use sort to ensure correct names used for predictors


######### maxent models ###############
spp_dirs <- list.dirs(model_dir, full.names=T, recursive=F)


######### thresholded predictions for each species #############

predictBinary <- function(model, predictors, threshold){
        # function to determine threshold, make model prediction, and return threshold prediction
        eval <- evaluate(model@presence, model@absence, model)
        thresh <- threshold(eval, stat=threshold)
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
                           p <- stack(lapply(files,function(x) readRDS(paste(clim_dir,x,sep="/"))))
                           names(p) <- climnames
                           bp <- predictBinary(m, p, maxent_threshold_stat)
                           saveRDS(bp, paste0(spp, "/BinaryRangePrediction.rds"))
                           return("success")
                   }

stopCluster(cl)
table(unlist(results))

######## richness map summed across species ###########

# exclude species that occur in less than our threshold for minimum number of cells
freqs <- read.csv(paste0(project_stem_dir, "/git_files/data/species_occurrence_counts.csv"), stringsAsFactors=F)
min_cells <- richness_min_cells
valid_species <- freqs$spp[freqs$ncells >= min_cells]

sumRasters <- function(dir, pattern, species){
        # function returns the sum of all raster files in 'dir' (and subdirs) that match 'pattern' and are among 'species'
        files <- list.files(dir, pattern=pattern, full.names=T, recursive=T)
        files <- files[basename(dirname(files)) %in% species]
        for(f in files[1:length(files)]){
                r <- readRDS(f)
                if(class(r)!="RasterLayer") next()
                if(f == files[1]) richness <- r else(richness <- richness + r)
        }
        return(richness)
}

richness <- sumRasters(model_dir, "BinaryRangePrediction.rds", valid_species)
writeRaster(richness, paste0(richness_dir, "/richness_810m_min", min_cells, ".tif"), overwrite=T)


######### plot #########

library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)

richness <- raster(paste0(richness_dir, "/richness_810m_min", min_cells, ".tif"))
r <- as.data.frame(rasterToPoints(richness))
names(r) <- c("x", "y", "layer")

# map
p <- ggplot(r, aes(x, y, fill=layer)) +
        geom_raster() +
        scale_fill_viridis() +
        coord_fixed(ratio=1) +
        theme(panel.background=element_blank(), panel.grid=element_blank(),
              axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
              legend.position="none", title=element_text(size=25)) +
        guides(fill = guide_colorbar(barwidth=20)) +
        labs(title=paste0("PLANT SPECIES RICHNESS\n(MaxEnt models for ", length(valid_species), 
                          " species with records in >= ", min_cells, " cells)"))

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

png(paste0(richness_dir, "/richness_min", min_cells, ".png"), width = 1000, height = 1500)
plot(p)
print(l, 
      vp=viewport(x = .5, y = .55, 
                  width = unit(0.45, "npc"), height = unit(0.3, "npc"),
                  just = c("left", "bottom")))
dev.off()



############## comparison of thresholds 0, 10, 30 ####################

r <- list.files(richness_dir, pattern="richness_810m", full.names=T)
r <- stack(r)
png(paste0(richness_dir, "/richness_thresholds_scatterplot.png"), width = 800, height = 800)
pairs(r)
dev.off()

v <- as.data.frame(rasterToPoints(r))
v[,3:5] <- apply(v[,3:5], 2, scale)
v$difference <- v$richness_810m_min30 - v$richness_810m_min0

p <- ggplot(v, aes(x, y, fill=difference)) +
        geom_raster() +
        scale_fill_gradient2(min="darkred", mid="gray", high="darkblue", midpoint=0) +
        coord_fixed(ratio=1) +
        theme(panel.background=element_blank(), panel.grid=element_blank(),
              axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
              legend.position="none", title=element_text(size=25)) +
        guides(fill = guide_colorbar(barwidth=20)) +
        labs(title="difference in standard scores, min30-min0")
h <- v %>%
        mutate(lyr = plyr::round_any(difference, .1)) %>%
        group_by(lyr) %>%
        summarize(n=n())
l <- ggplot(h, aes(lyr, n, fill=lyr)) + 
        geom_bar(stat="identity", width=.1) + 
        scale_fill_gradient2(min="darkred", mid="gray", high="darkblue", midpoint=0) +
        theme(panel.background=element_blank(), panel.grid=element_blank(),
              axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none") +
        labs(x="z(min30) - z(min0)")

png(paste0(richness_dir, "/richness_difference_30_0.png"), width = 1000, height = 1500)
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

cl <- makeCluster(nodes)
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

writeRaster(sumRasters(model_dir, "BinaryRangePrediction_25k.rds", valid_species), 
            paste0(richness_dir, "/richness_25k_min", min_cells, ".tif"), overwrite=T)
writeRaster(sumRasters(model_dir, "BinaryRangePrediction_50k.rds", valid_species),
            paste0(richness_dir, "/richness_50k_min", min_cells, ".tif"), overwrite=T)


### compare to occurrence density maps ###

for(res in c("25", "50")){
        r <- stack(raster(paste0(richness_dir, "/richness_min", min_cells, "_", res, "k.tif")),
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
