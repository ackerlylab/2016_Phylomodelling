

# Threshold Maxent predictions and sum to richness
# Matthew Kling
# January 2016


library(raster)
library(dismo)
library(doParallel)
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)


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

for(min_cells in c(0,10,30)){ # calculte richness for each of the three thresholds
        
        # exclude species that occur in less than our threshold for minimum number of cells
        freqs <- read.csv(paste0(project_stem_dir, "/git_files/data/species_occurrence_counts.csv"), stringsAsFactors=F)
        valid_species <- freqs$spp[freqs$ncells >= min_cells]
        
        # compute richness
        richness <- sumRasters(model_dir, "BinaryRangePrediction.rds", valid_species)
        
        # clip to state border
        cali <- rgdal::readOGR("C:/Lab_projects/2016_Phylomodelling/Data/Shapefiles/states", "cb_2014_us_state_500k")
        cali <- cali[cali$NAME=="California",]
        cali <- spTransform(cali, crs(richness))
        richness <- mask(richness, cali)
        
        # save raster file
        writeRaster(richness, paste0(richness_dir, "/richness_810m_min", min_cells, ".tif"), overwrite=T)
}

######### plots #########



for(min_cells in c(0,10,30)){ # create richness maps for each of the three thresholds
        
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
        
        png(paste0(richness_dir, "/richness_min", min_cells, ".png"), width = 1200, height = 1500)
        plot(p)
        print(l, 
              vp=viewport(x = .5, y = .6, 
                          width = unit(0.4, "npc"), height = unit(0.25, "npc"),
                          just = c("left", "bottom")))
        dev.off()
}




################# upscale to 25 and 50 km ###############

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
for(min_cells in c(0,10,30)){
        freqs <- read.csv(paste0(project_stem_dir, "/git_files/data/species_occurrence_counts.csv"), stringsAsFactors=F)
        valid_species <- freqs$spp[freqs$ncells >= min_cells]
        
        writeRaster(sumRasters(model_dir, "BinaryRangePrediction_25k.rds", valid_species),
                    paste0(richness_dir, "/richness_25k_min", min_cells, ".tif"), overwrite=T)
        writeRaster(sumRasters(model_dir, "BinaryRangePrediction_50k.rds", valid_species),
                    paste0(richness_dir, "/richness_50k_min", min_cells, ".tif"), overwrite=T)
}




############## comparison of THRESHOLDS: 0, 10, 30 ####################

### maxent layers ###
for(res in c("25k", "50k", "810m")){
        
        r <- list.files(richness_dir, pattern=paste0("richness_", res), full.names=T)
        r <- stack(r)
        r <- mask(r, cali)
        
        png(paste0(richness_dir, "/richness_thresholds_scatterplot_", res, "_maxent.png"), width = 800, height = 800)
        pairs(r)
        dev.off()
        
        v <- as.data.frame(rasterToPoints(r))
        v[,3:5] <- apply(v[,3:5], 2, scale)
        v$difference <- v[,5] - v[,3]
        
        p <- ggplot(v, aes(x, y, fill=difference)) +
                geom_raster() +
                scale_fill_gradient2(min="darkred", mid="gray", high="darkblue", midpoint=0) +
                coord_fixed(ratio=1) +
                theme(panel.background=element_blank(), panel.grid=element_blank(),
                      axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                      legend.position="none", title=element_text(size=25)) +
                guides(fill = guide_colorbar(barwidth=20)) +
                labs(title=paste0("Difference in relative richness,\nno minimum vs 30-cell minimum,\n", res, " maxent"))
        binwidth <- range(v$difference, na.rm=T)
        binwidth <- (binwidth[2]-binwidth[1])/20
        h <- v %>%
                mutate(lyr = plyr::round_any(difference, binwidth)) %>%
                group_by(lyr) %>%
                summarize(n=n())
        l <- ggplot(h, aes(lyr, n, fill=lyr)) + 
                geom_bar(stat="identity", width=binwidth) + 
                scale_fill_gradient2(min="darkred", mid="gray", high="darkblue", midpoint=0) +
                theme(panel.background=element_blank(), panel.grid=element_blank(),
                      axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none") +
                labs(x="z(min30) - z(min0)")
        
        png(paste0(richness_dir, "/richness_difference_30_0_", res, "_maxent.png"), width = 1200, height = 1500)
        plot(p)
        print(l, 
              vp=viewport(x = .5, y = .6, 
                          width = unit(0.4, "npc"), height = unit(0.25, "npc"),
                          just = c("left", "bottom")))
        dev.off()
}


### and point richness layers ###
for(res in c("25", "50")){
        for(stat in c("raw", "Rarefied")){
                
                r <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Richness/point_richness_rasters", pattern=stat, full.names=T)
                r <- r[grepl(".grd", r) & grepl(paste0(res, "k"), r)]
                layers <- basename(r)
                r <- stack(r)
                names(r) <- layers
                
                # pairs plot
                png(paste0(richness_dir, "/richness_thresholds_scatterplot_", res, "k_", stat, ".png"), width = 800, height = 800)
                pairs(r)
                dev.off()
                
                # difference in relative richness, 30 vs 0 thresholds
                v <- as.data.frame(rasterToPoints(r))
                v[,3:6] <- apply(v[,3:6], 2, scale)
                v$difference <- v[,6] - v[,4]
                p <- ggplot(v, aes(x, y, fill=difference)) +
                        geom_raster() +
                        scale_fill_gradient2(min="darkred", mid="gray", high="darkblue", midpoint=0) +
                        coord_fixed(ratio=1) +
                        theme(panel.background=element_blank(), panel.grid=element_blank(),
                              axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                              legend.position="none", title=element_text(size=25)) +
                        guides(fill = guide_colorbar(barwidth=20)) +
                        labs(title=paste0("difference in relative richness,\nno minimum vs 30-cell minimum,\n", res, "k ", stat))
                binwidth <- range(v$difference, na.rm=T)
                binwidth <- (binwidth[2]-binwidth[1])/20
                h <- v %>%
                        mutate(lyr = plyr::round_any(difference, binwidth)) %>%
                        group_by(lyr) %>%
                        summarize(n=n())
                l <- ggplot(h, aes(lyr, n, fill=lyr)) + 
                        geom_bar(stat="identity", width=binwidth) + 
                        scale_fill_gradient2(min="darkred", mid="gray", high="darkblue", midpoint=0) +
                        theme(panel.background=element_blank(), panel.grid=element_blank(),
                              axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none") +
                        labs(x="z(min30) - z(min0)")
                
                png(paste0(richness_dir, "/richness_difference_30_0_", res, "k_", stat, ".png"), width = 1200, height = 1500)
                plot(p)
                print(l, 
                      vp=viewport(x = .5, y = .6, 
                                  width = unit(0.4, "npc"), height = unit(0.25, "npc"),
                                  just = c("left", "bottom")))
                dev.off()
                
        }
}




############## comparison of STATISTICS: maxent, raw, rarefied ####################

resolutions <- c(25, 50)
thresholds <- c(0, 10, 30)

for(res in resolutions){
        for(thresh in thresholds){
                
                # load and stack rasters for the three stats
                thresh_tag <- paste0("n", thresh); if(thresh==0) thresh_tag <- "all_bcm"
                r <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Richness/point_richness_rasters", 
                                pattern=thresh_tag, full.names=T)
                r <- r[grepl(".grd", r) & grepl(paste0(res, "k"), r)]
                m <- list.files(richness_dir, pattern=".tif", full.names=T)
                m <- m[grepl(paste0(res, "k"), m) & grepl(paste0("min", thresh), m)]
                s <- stack(c(r, m))
                names(s) <- c("rarefied", "raw", "maxent")
                s <- mask(s, s$raw) # mask maxent to CA domain
                
                # pairs plot
                png(paste0(richness_dir, "/richness_statistics_scatterplot_", res, "k_min", thresh, ".png"), width = 800, height = 800)
                pairs(s)
                dev.off()
                
                # basic richness maps
                s <- as.data.frame(rasterToPoints(s)) %>%
                        gather(stat, value, -x, -y) %>%
                        mutate(stat=factor(stat, levels=c("raw", "rarefied", "maxent")))
                p1 <- ggplot(s, aes(x, y, fill=value)) + 
                        geom_raster() + 
                        facet_grid(.~stat) +
                        scale_fill_viridis(na.value="white")  +
                        coord_fixed(ratio=1) +
                        theme(panel.background=element_blank(), panel.grid=element_blank(),
                              axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                              legend.position="top", title=element_text(size=25), strip.text=element_text(size=25)) +
                        guides(fill = guide_colorbar(barwidth=20)) +
                        labs(fill="richness (# spp)  ")
                
                # normalized richness maps
                t <- s %>%
                        group_by(stat) %>%
                        mutate(value = scale(value))
                #mutate(value = value/max(value, na.rm=T))
                p2 <- ggplot(t, aes(x, y, fill=value)) + 
                        geom_raster() + 
                        facet_grid(.~stat) +
                        scale_fill_viridis(na.value="white")  +
                        coord_fixed(ratio=1) +
                        theme(panel.background=element_blank(), panel.grid=element_blank(),
                              axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                              legend.position="top", title=element_text(size=25), strip.text=element_text(size=25)) +
                        guides(fill = guide_colorbar(barwidth=20)) +
                        labs(fill="normalized richness (z-score)  ")
                
                # differnce maps from maxent
                u <- t %>%
                        spread(stat, value) %>%
                        mutate(rarefied=rarefied-maxent, raw=raw-maxent, maxent=maxent-maxent) %>%
                        gather(stat, value, -x, -y)
                mag <- max(abs(range(u$value, na.rm=T)))
                p3 <- ggplot(u, aes(x, y, fill=value)) + 
                        geom_raster() + 
                        facet_grid(.~stat) +
                        scale_fill_gradientn(colours=c("darkblue", "dodgerblue", "gray", "orange", "darkred"), 
                                             na.value="white", limits=c(-mag, mag))  +
                        coord_fixed(ratio=1) +
                        theme(panel.background=element_blank(), panel.grid=element_blank(),
                              axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                              legend.position="top", title=element_text(size=25), strip.text=element_text(size=25)) +
                        guides(fill = guide_colorbar(barwidth=20)) +
                        labs(fill="sample minus maxent (normalized)  ")
                
                # save combo plot
                p <- arrangeGrob(textGrob(label=paste0("Species richness: statistics compared (", res, "km, ", "min ", thresh, " cells)"),
                                          gp=gpar(fontsize=30)), 
                                 p1, p2, p3, ncol=1, heights=c(.2, 1, 1, 1))
                png(paste0(richness_dir, "/richness_comparison_", res, "k_min", thresh, ".png"), width=1000, height=1500)
                grid.draw(p)
                dev.off()
                
                # difference from maxent as function of raw
                z <- s %>%
                        spread(stat, value) %>%
                        gather(stat, value, rarefied, raw) %>%
                        group_by(stat) %>%
                        mutate(maxent_scaled=scale(maxent), value_scaled=scale(value),
                               maxent_pmax=maxent/max(na.omit(maxent)), value_pmax=value/max(na.omit(value)),
                               diff_basic=value-maxent,
                               diff_scaled=value_scaled-maxent_scaled,
                               diff_pmax=value_pmax-maxent_pmax) %>%
                        select(x, y, stat, value, diff_basic:diff_pmax) %>%
                        gather(diff_stat, diff_value, diff_basic:diff_pmax)
                        
                plt <- ggplot(z, aes(value, diff_value)) + 
                        geom_point() +
                        geom_smooth(se=F) +
                        facet_grid(diff_stat~stat, scales="free") +
                        theme_bw() +
                        theme(legend.position="top") +
                        labs(x="richness", y="",
                             title=paste0("sample richness minus maxent richness, ", res, "km, ", "min ", thresh, " cells"))
                ggsave(paste0(richness_dir, "/error_saturation_scatterplot_", res, "k_min", thresh, ".png"), plt, width=8, height=10, units="in")
        }
}




######### add range sizes to species table ############

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
names(rd) <- c("spp", "maxent810m", "maxent25km", "maxent50km")
rd <- full_join(freqs, rd)
rd <- select(rd, -X)
rd <- select(rd, spp, nrecords, ncells, maxent810m, maxent25km, maxent50km)
write.csv(rd, paste0(project_stem_dir, "/git_files/data/species_occurrence_counts.csv"), row.names=F)

# pairs plots, linear and log scales
png(paste0(richness_dir, "/species_range_scatterplot.png"), width=1000, height=1000)
pairs(rd[,2:6], cex=.1)
dev.off()
png(paste0(richness_dir, "/species_range_scatterplot_loglog.png"), width=1000, height=1000)
pairs(log10(rd[,2:6]), cex=.1)
dev.off()
