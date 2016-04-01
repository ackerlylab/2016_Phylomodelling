

## to do
# export raster layers separately
# add unfiltered cch records as additional panel



library(raster)
library(dismo)
library(doParallel)
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)
library(rgdal)


clim_dir <- filled_climate_data_dir
model_dir <- paste0(project_stem_dir, "/Output/Maxent/v5")
richness_dir <- paste0(project_stem_dir, "/Output/Richness/V5")
pr_dir <- paste0(project_stem_dir, "/Data/Richness/point_richness_rasters")
spdir <- occurrence_data_dir_processed
maxent_background <- paste(spdir, 'combined/10000_CA_plants_bg_810m_occbias.rdata', sep='')


outdir <- "C:/Lab_projects/2016_Phylomodelling/Output/Maxent/range_constraint_test"

# load occurrences
allocc <- readRDS("C:/Lab_projects/2016_Phylomodelling/Data/Species/occurrences_clean.rds")
species <- unique(allocc$current_name_binomial)

# load ecoregion shapefile
ecoregions <- readOGR("C:/Lab_projects/2016_Phylomodelling/Data/Shapefiles/jepson_ecoregions", "jepson_ecoregions")
proj4string(ecoregions) <- "+proj=longlat +ellps=WGS84"
ecoregions <- spTransform(ecoregions, crs(allocc))


# climate data
files <- list.files(path=clim_dir, pattern='810m_filled3x', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt"))
p <- stack(lapply(files,function(x) readRDS(paste(clim_dir,x,sep="/"))))
names(p) <- climnames

# california boundary
cali <- rgdal::readOGR("C:/Lab_projects/2016_Phylomodelling/Data/Shapefiles/states", "cb_2014_us_state_500k")
cali <- cali[cali$NAME=="California",]
cali <- spTransform(cali, crs(readRDS(paste0(clim_dir, "/", files[1]))))



spp_dirs <- list.dirs(model_dir, full.names=T, recursive=F)




binary_range <- function(points=NULL, background=NULL, maxmod=NULL, predictors=NULL, threshold=NULL, 
                         constraint=NULL, constraint_method="none", sigma=50, continuous=FALSE){
        
        pred <- predict(maxmod, predictors)
        if(continuous) return(pred)
        
        calrst <- predictors[[1]]
        values(calrst)[!is.na(values(calrst))] <- 0
        
        
        if(constraint_method=="distance"){
                if(class(constraint)=="SpatialPolygonsDataFrame"){
                        dst <- mask(calrst, constraint)
                        dst <- distance(dst) # this is the slow step
                } else if(class(constraint)=="SpatialPointsDataFrame"){
                        dst <- distanceFromPoints(calrst, constraint)
                } else{
                        stop("invalid constraint feature for distance method; must be spatial points or polygons")
                }
                
                dst <- mask(dst, calrst)
                gauss <- function(x, sigma) exp(-.5 * (x/sigma) ^ 2)
                dst <- calc(dst, function(x) gauss(x, sigma*1000))
                pred <- pred * dst
        }
        
        extract <- raster::extract
        eval <- evaluate(p=extract(pred, points), a=extract(pred, background))
        thresh <- threshold(eval, stat=threshold)
        pred <- reclassify(pred, c(0, thresh, 0, thresh, 1, 1))
        
        if(constraint_method=="clip"){
                pred <- mask(pred, constraint)
                pred <- mask(sum(pred, calrst, na.rm=T), calrst)
        }
        
        return(pred)
}




point_buffs <- list.files(paste0(project_stem_dir, "/Output/Range_polygons/occurrence_buffers"), full.names=T)
hull_buffs <- list.files(paste0(project_stem_dir, "/Output/Range_polygons/buffered_hulls"), full.names=T)
predictors <- p
background <- readRDS(maxent_background)
threshold <- maxent_threshold_stat
regions <- ecoregions
sigma <- 50
cali_mask <- mask(p[[1]], cali)



# raw (pre-cleaning) occurrences
get_raw_occ <- function(x){
        file <- paste0(x, "/Suitability_data.rdata")
        if(file.exists(file)){
                x <- readRDS(paste0(x, "/Suitability_data.rdata"))
                as.data.frame(x$occur)
        }
}
raw_occ_dirs <- list.dirs(paste0(maxent_output_dir, "/V2"), recursive=F)




##############################################
#### v2: modified after consult w david
##############################################

myspp <- spp_dirs

cl <- makeCluster(nodes)
registerDoParallel(cl)
results <- foreach(spp = myspp,
                   .packages=c("raster", "dismo", "ggplot2", "tidyr", 
                               "dplyr", "grid", "gridExtra", "rgdal")) %dopar% {
                                       
                                       ###### load species-specific data #####
                                       
                                       points <- allocc[allocc$current_name_binomial==basename(spp),]
                                       maxmod <- readRDS(paste0(spp, "/ModelObject.rdata"))
                                       if(class(maxmod)=="try-error") return("aborted -- no maxent model")
                                       
                                       
                                       ##### generate ranges ######
                                       
                                       # continuous maxent
                                       cmx <- binary_range(maxmod=maxmod, predictors=predictors, continuous=T)
                                       
                                       # straight maxent
                                       mxm <- binary_range(points, background, maxmod, predictors, threshold)
                                      
                                       # clip to ecoregion
                                       ecoregions <- over(regions, points)
                                       if(all(is.na(ecoregions$id))) return("aborted -- no records within ecoregions")
                                       ecoregions <- regions[!is.na(ecoregions$id),]
                                       erc <- binary_range(points, background, maxmod, predictors, threshold,
                                                           constraint=ecoregions, constraint_method="clip")
                                      
                                       # distance from points
                                       ptd <- binary_range(points, background, maxmod, predictors, threshold,
                                                           constraint=points, constraint_method="distance")
                                       
                                       # combine and save
                                       s <- stack(cmx, mxm, erc, ptd)
                                       names(s) <- c("continuous_maxent", "binary_maxent",
                                                     "ecoregion_clip", "distance_hybrid")
                                       s <- mask(s, cali_mask)
                                       writeRaster(s, paste0(outdir, "/v2/rasters/multiband/", basename(spp), ".tif"), overwrite=T)
                                       
                                       # save individual layers
                                       writeRaster(s[[1]], paste0(outdir, "/v2/rasters/continuous_maxent/", basename(spp), ".tif"), overwrite=T)
                                       writeRaster(s[[2]], paste0(outdir, "/v2/rasters/binary_maxent/", basename(spp), ".tif"), overwrite=T)
                                       writeRaster(s[[3]], paste0(outdir, "/v2/rasters/ecoregion_clip/", basename(spp), ".tif"), overwrite=T)
                                       writeRaster(s[[4]], paste0(outdir, "/v2/rasters/distance_hybrid/", basename(spp), ".tif"), overwrite=T)
                                       
                                       
                                       ##### prep plot data ######
                                       s <- aggregate(s, 2, fun=max, na.rm=T) # aggregate for plotting speed & memory
                                       d <- as.data.frame(rasterToPoints(s)) %>%
                                               gather(method, presence, -x, -y)
                                       
                                       pd <- d[d$method=="ecoregion_clip",]
                                       pd$method <- "clean_occurrences"
                                       pd$presence <- 0
                                       d <- rbind(d, pd)
                                       pd$method <- "all_occurrences"
                                       d <- rbind(d, pd)
                                       
                                       pd <- as.data.frame(coordinates(points))
                                       names(pd) <- c("x", "y")
                                       pd$method <- "clean_occurrences"
                                       
                                       raw_occ <- try(get_raw_occ(raw_occ_dirs[grepl(basename(spp), raw_occ_dirs)]) %>%
                                               select(longitude, latitude))
                                       if(class(raw_occ)=="try-error") return("aborted -- unclean occurrencse not found")
                                       raw_occ$method <- "all_occurrences"
                                       
                                       variables <- c("all_occurrences", "clean_occurrences", 
                                                      "continuous_maxent", "binary_maxent",
                                                      "ecoregion_clip", "distance_hybrid")
                                       
                                       d$method <- factor(d$method, levels=variables)
                                       pd$method <- factor(pd$method, levels=variables)
                                       raw_occ$method <- factor(raw_occ$method, levels=variables)
                                       
                                       d$continuous[d$method=="continuous_maxent"] <- d$presence[d$method=="continuous_maxent"]
                                       d$presence[d$method=="continuous_maxent"] <- NA
                                       
                                       best_occ <- as.data.frame(points[which.max(raster::extract(cmx, points))[1],])
                                       best_occ$method <- "clean_occurrences"
                                       best_occ$method <- factor(best_occ$method, levels=variables)
                                       
                                       best_mx <- filter(d, continuous==max(continuous, na.rm=T)) %>%
                                               sample_n(1)
                                       best_mx$method <- "binary_maxent"
                                       best_mx$method <- factor(best_mx$method, levels=variables)
                                       
                                       
                                       
                                       ##### build plots #####
                                       
                                       eb <- ggplot2::element_blank()
                                       plt <- ggplot() +
                                               geom_raster(data=d, aes(x, y), fill="gray80") +
                                               geom_raster(data=d[d$presence==1 & !is.na(d$presence),], aes(x, y), fill="darkred") +
                                               geom_raster(data=d[d$method=="continuous_maxent",], aes(x, y, fill=continuous)) +
                                               geom_point(data=pd, aes(x, y), color="darkred", size=3, shape=3) +
                                               geom_point(data=raw_occ, aes(longitude, latitude), color="darkred", size=3, shape=3) +
                                               geom_point(data=best_occ, aes(longitude, latitude), color="yellow", size=3, shape=10) +
                                               geom_point(data=best_mx, aes(x, y), color="yellow", size=3, shape=10) +
                                               scale_fill_gradientn(colours=c("gray80", "darkred")) +
                                               facet_wrap(~method, nrow=1, labeller=as_labeller(function(x) paste0(sub("_", "\n ", x), "\n"))) +
                                               coord_fixed(ratio=1.3, expand=c(0,0)) +
                                               labs(title=paste0(basename(spp), "\n")) +
                                               theme(panel.background=eb, panel.grid=eb, 
                                                     axis.text=eb, axis.ticks=eb, axis.title=eb,
                                                     strip.background=eb, legend.position="none", 
                                                     title=element_text(color="darkred", size=30), 
                                                     strip.text=element_text(color="darkred", size=14),
                                                     panel.margin = unit(0, "lines"))
                                       
                                       ggsave(paste0(outdir, "/v2/charts/", basename(spp), ".png"), plt, width=12, height=9)
                               }
stopCluster(cl)


stop("finito")

##############################################
#### v1
##############################################
myspp <- sample(spp_dirs, 1000)

cl <- makeCluster(nodes)
registerDoParallel(cl)
results <- foreach(spp = myspp,
                   .packages=c("raster", "dismo", "ggplot2", "tidyr", 
                               "dplyr", "grid", "gridExtra", "rgdal")) %dopar% {
                                       
                                       
                                       if(file.exists(paste0(outdir, "/rasters/", basename(spp)))) return("skipped -- outfile already exists")
                                       
                                       ###### load species-specific data #####
                                       
                                       points <- allocc[allocc$current_name_binomial==basename(spp),]
                                       maxmod <- readRDS(paste0(spp, "/ModelObject.rdata"))
                                       if(class(maxmod)=="try-error") return("aborted -- no maxent model")
                                       
                                       
                                       ##### generate ranges ######
                                       
                                       # straight maxent
                                       mxm <- binary_range(points, background, maxmod, predictors, threshold)
                                       
                                       # clip to point buffer
                                       buff <- try(readRDS(point_buffs[sub("\\.rds", "", basename(point_buffs)) == basename(spp)]))
                                       if(class(buff)=="try-error") return("aborted -- no point buffer shapefile")
                                       buff <- spTransform(buff, crs(predictors))
                                       pbc <- binary_range(points, background, maxmod, predictors, threshold,
                                                           constraint=buff, constraint_method="clip")
                                       
                                       # clip to convex hull buffer
                                       buff <- try(readRDS(hull_buffs[sub("\\.rds", "", basename(hull_buffs)) == basename(spp)]))
                                       if(class(buff)=="try-error") return("aborted -- no hull buffer shapefile")
                                       buff <- spTransform(buff, crs(predictors))
                                       hbc <- binary_range(points, background, maxmod, predictors, threshold,
                                                           constraint=buff, constraint_method="clip")
                                       
                                       # clip to ecoregion
                                       ecoregions <- over(regions, points)
                                       ecoregions <- regions[!is.na(ecoregions$id),]
                                       erc <- binary_range(points, background, maxmod, predictors, threshold,
                                                           constraint=ecoregions, constraint_method="clip")
                                       
                                       # distance from ecoregion
                                       erd <- binary_range(points, background, maxmod, predictors, threshold,
                                                           constraint=ecoregions, constraint_method="distance")
                                       
                                       # distance from points
                                       ptd <- binary_range(points, background, maxmod, predictors, threshold,
                                                           constraint=points, constraint_method="distance")
                                       
                                       # combine and save
                                       s <- stack(mxm, pbc, hbc, erc, erd, ptd)
                                       names(s) <- c("straight_maxent", "point_buffer_clip", "hull_buffer_clip",
                                                     "ecoregion_clip", "ecoregion_distance", "points_distance")
                                       saveRDS(s, paste0(outdir, "/rasters/multiband/", basename(spp)))
                                       
                                       
                                       
                                       ##### make plots ######
                                       s <- aggregate(s, 2, fun=max, na.rm=T)
                                       d <- as.data.frame(rasterToPoints(s)) %>%
                                               gather(method, presence, -x, -y)
                                       
                                       pd <- d[d$method=="ecoregion_clip",]
                                       pd$method <- "occurrences"
                                       pd$presence <- 0
                                       d <- rbind(d, pd)
                                       
                                       pd <- as.data.frame(coordinates(points))
                                       names(pd) <- c("x", "y")
                                       pd$method <- "occurrences"
                                       
                                       d$method <- factor(d$method, levels=c("occurrences", "straight_maxent", "point_buffer_clip", "hull_buffer_clip",
                                                                             "ecoregion_clip", "ecoregion_distance", "points_distance"))
                                       pd$method <- factor(pd$method, levels=c("occurrences", "straight_maxent", "point_buffer_clip", "hull_buffer_clip",
                                                                               "ecoregion_clip", "ecoregion_distance", "points_distance"))
                                       
                                       eb <- ggplot2::element_blank()
                                       plt <- ggplot() +
                                               geom_raster(data=d, aes(x, y, fill=factor(presence))) +
                                               geom_point(data=pd, aes(x, y), color="darkred", size=3, shape=3) +
                                               scale_fill_manual(values=c("gray80", "darkred")) +
                                               facet_wrap(~method, nrow=2) +
                                               coord_fixed(ratio=1.3) +
                                               labs(title=paste0(basename(spp), "\n")) +
                                               theme(panel.background=eb, panel.grid=eb, 
                                                     axis.text=eb, axis.ticks=eb, axis.title=eb,
                                                     strip.background=eb, legend.position="none", 
                                                     title=element_text(color="darkred", size=30), 
                                                     strip.text=element_text(color="darkred", size=14))
                                       
                                       ggsave(paste0(outdir, "/charts/", basename(spp), ".png"), plt, width=12, height=9)
                               }
stopCluster(cl)


