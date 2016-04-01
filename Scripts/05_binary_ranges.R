

# Threshold Maxent predictions and upscale to coarser resolutions
# Matthew Kling
# Last updated March 2016


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
model_dir <- paste0(project_stem_dir, "/Output/Maxent/v5")
richness_dir <- paste0(project_stem_dir, "/Output/Richness/V5")
pr_dir <- paste0(project_stem_dir, "/Data/Richness/point_richness_rasters")


######### climate data ###########
files <- list.files(path=clim_dir, pattern='810m_filled3x', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt")) #use sort to ensure correct names used for predictors


######### maxent models ###############
spp_dirs <- list.dirs(model_dir, full.names=T, recursive=F)

######### california boundary ############
cali <- rgdal::readOGR("C:/Lab_projects/2016_Phylomodelling/Data/Shapefiles/states", "cb_2014_us_state_500k")
cali <- cali[cali$NAME=="California",]
cali <- spTransform(cali, crs(readRDS(paste0(clim_dir, "/", files[1]))))


######### thresholded predictions for each species #############

predictBinary <- function(model, predictors, threshold){
        # function to determine threshold, make model prediction, and return threshold prediction
        eval <- evaluate(model@presence, model@absence, model)
        thresh <- threshold(eval, stat=threshold)
        pred <- predict(model, predictors)
        pred <- reclassify(pred, c(0, thresh, 0, thresh, 1, 1))
        return(pred)
}

# occurrence buffer shapefiles
#buffs <- list.files(paste0(project_stem_dir, "/Output/Range_polygons/occurrence_buffers"), full.names=T)
buffs <- list.files(paste0(project_stem_dir, "/Output/Range_polygons/buffered_hulls"), full.names=T)

cl <- makeCluster(nodes)
registerDoParallel(cl)
results <- foreach(spp = spp_dirs,
                   .packages=c("raster", "dismo")) %dopar% {
                           
                           m <- readRDS(paste0(spp, "/ModelObject.rdata"))
                           if(class(m)=="try-error") return("no maxent model")
                           p <- stack(lapply(files,function(x) readRDS(paste(clim_dir,x,sep="/"))))
                           names(p) <- climnames
                           bp <- predictBinary(m, p, maxent_threshold_stat)
                           
                           # clip to CA boundary and save
                           bp <- mask(bp, cali)
                           saveRDS(bp, paste0(spp, "/BinaryRangePrediction.rds"))
                           
                           # save a version clipped to point buffer
                           buff <- buffs[sub("\\.rds", "", basename(buffs)) == basename(spp)]
                           if(length(buff)==0) return("no buffer shapefile")
                           buff <- readRDS(buff)
                           buff <- spTransform(buff, crs(bp))
                           bp <- mask(bp, buff)
                           #saveRDS(bp, paste0(spp, "/BufferClippedMaxent.rds"))
                           saveRDS(bp, paste0(spp, "/HullBufferClippedMaxent.rds"))
                           
                           return("success")
                   }
stopCluster(cl)
table(unlist(results))




################# upscale to 25 and 50 km ###############

upscale <- function(file, template, tag){
        outfile <- sub("\\.rds", paste0("_", tag, ".rds"), file)
        #if(file.exists(outfile)) return("file already exists")
        f <- readRDS(file)
        f <- as.data.frame(rasterToPoints(f))
        coordinates(f) <- c("x", "y")
        f <- rasterize(f, template, field="layer", fun=function(x, ...){max(na.omit(x))})
        saveRDS(f, outfile)
        return("upscaled")
}

pr_files <- list.files(pr_dir, pattern=".grd", full.names=T)

for(rangetype in c("BinaryRangePrediction.rds", "BufferClippedMaxent.rds")){
        cl <- makeCluster(nodes)
        registerDoParallel(cl)
        results <- foreach(spp = paste0(spp_dirs, "/", rangetype),
                           .packages=c("raster")) %dopar% {
                                   if(!file.exists(spp)) return("no maxent model")
                                   upscale(spp, template=raster(pr_files[1]), "25k")
                                   upscale(spp, template=raster(pr_files[2]), "50k")
                                   return("success")
                           }
        stopCluster(cl)
}






######## export png range maps for each species ###########

# get correct projection
#prj <- crs(readRDS(paste0(spp_dirs[1], "/BinaryRangePrediction.rds")))

# load occurrences
allocc <- readRDS("C:/Lab_projects/2016_Phylomodelling/Data/Species/occurrences_clean.rds")

# and buffer polygons
#buffs <- list.files("C:/Lab_projects/2016_Phylomodelling/Output/Range_polygons/occurrence_buffers", full.names=T)
buffs <- list.files("C:/Lab_projects/2016_Phylomodelling/Output/Range_polygons/buffered_hulls", full.names=T)

# geom_holygon function from http://qiita.com/kohske/items/9272e29a75d32416ff5e fixes polygon hole bug in ggplot2::geom_polygon
library(ggplot2)
library(proto)
library(grid)
GeomHolygon <- ggproto(
        "GeomHolygon", 
        GeomPolygon,
        extra_params = c("na.rm", "rule"),
        draw_panel = function(data, scales, coordinates, rule) {
                n <- nrow(data)
                if (n == 1) 
                        return(zeroGrob())
                
                munched <- coord_munch(coordinates, data, scales)
                munched <- munched[order(munched$group), ]
                
                first_idx <- !duplicated(munched$group)
                first_rows <- munched[first_idx, ]
                
                ggplot2:::ggname(
                        "geom_holygon", 
                        pathGrob(munched$x, munched$y, default.units = "native", 
                                 id = munched$group, rule = rule, 
                                 gp = gpar(col = first_rows$colour, 
                                           fill = alpha(first_rows$fill, first_rows$alpha), 
                                           lwd = first_rows$size * .pt, 
                                           lty = first_rows$linetype)))
        }
)
geom_holygon <- function (mapping = NULL, data = NULL, stat = "identity", position = "identity", 
                          na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, rule = "winding", ...) {
        ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomHolygon, 
                       position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
                       params = list(na.rm = na.rm , rule = rule, ...))
}

# function to generate chart
rangemap <- function(dir){
        files <- list.files(dir, pattern="BinaryRangePrediction", full.names=T)
        if(length(files)==0) return("no maxent prediction")
        rds2df <- function(x) as.data.frame(rasterToPoints(readRDS(x)))
        r <- lapply(files, rds2df)
        for(i in 1:length(r)) r[[i]]$resolution <- c("810m", "25km", "50km")[i]
        r <- do.call("rbind", r)
        names(r)[3] <- "presence"
        
        # add point occurrences
        occ <- allocc[allocc$current_name_binomial==basename(dir),]
        occ <- as.data.frame(occ)
        occ$resolution <- "records"
        occ$resolution_f <- factor(occ$resolution, levels=c("records", "810m", "25km", "50km"))
        occbg <- r[r$resolution=="810m",]
        occbg$resolution <- "records"
        occbg$presence <- 0
        r <- rbind(r, occbg)
        r$resolution_f <- factor(r$resolution, levels=c("records", "810m", "25km", "50km"))
        
        # add occurrence buffer
        buff <- buffs[sub("\\.rds", "", basename(buffs))==basename(dir)]
        buff <- readRDS(buff)
        buff <- spTransform(buff, crs(cali))
        buff <- gIntersection(buff, cali)
        buff <- fortify(buff)
        buff$resolution <- "810m"
        buff$resolution_f <- factor(buff$resolution, levels=c("records", "810m", "25km", "50km"))
        
        eb <- ggplot2::element_blank()
        p <- ggplot() +
                geom_raster(data=r, aes(x, y, fill=factor(presence))) +
                geom_point(data=occ, aes(longitude, latitude), 
                           color="darkred", shape=3, size=3) +
                geom_holygon(data=buff, aes(long, lat, group=group, order=order), 
                             color="orange", fill=NA, size=1) +
                scale_fill_manual(values=c("gray80", "darkred")) +
                facet_wrap(~resolution_f, nrow=1) +
                coord_fixed(ratio=1.3) +
                labs(title=paste0(basename(dir), "\n")) +
                theme(panel.background=eb, panel.grid=eb, 
                      axis.text=eb, axis.ticks=eb, axis.title=eb,
                      strip.background=eb, legend.position="none", 
                      title=element_text(color="darkred", size=30), 
                      strip.text=element_text(color="gray70", size=20))
        ggsave(paste0(dir, "/", basename(dir), " maxent rangemap.png"), 
               p, width=12, height=9, units="in")
}

cl <- makeCluster(nodes)
registerDoParallel(cl)
results <- foreach(spp = spp_dirs,
                   .packages=c("raster", "ggplot2", "rgeos", "sp", "grid", "proto")) %dopar% {
                           rangemap(spp)
                   }
stopCluster(cl)

# consolidate copies of these maps in a single location
maps <- list.files(dirname(spp_dirs[1]), recursive=T, full.names=T, pattern="maxent rangemap")
newmaps <- paste0("C:/Lab_projects/2016_Phylomodelling/Output/Charts/rangemaps_V5/", basename(maps))
file.copy(maps, newmaps, overwrite=T)


