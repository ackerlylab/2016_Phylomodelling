

library(doParallel)
library(rgeos)
library(dplyr)
library(raster)


cdir <- filled_climate_data_dir
outdir <- paste0(project_stem_dir, "/Output/Range_polygons")


# get CRS
proj <- crs(raster("C:/Lab_projects/2016_Phylomodelling/Output/Richness/V5/richness_810m_min0.tif"))

# load occurrences
allocc <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Species/processed3/atomic", full.names=T)
allocc <- lapply(allocc, readRDS)
allocc <- do.call("rbind", allocc)
allocc <- allocc[!is.na(allocc$longitude + allocc$latitude),]
coordinates(allocc) <- c("longitude", "latitude")
projection(allocc) <- '+proj=longlat +ellps=WGS84'
allocc <- spTransform(allocc, proj)

# exclude occurrences with no climate data
clim <- readRDS(list.files(cdir, pattern="cwd_", full.names=T))
allocc <- allocc[!is.na(raster::extract(clim, allocc)),]
species <- unique(allocc$current_name_binomial)

# load california boundary, to intersect with hulls
cali <- rgdal::readOGR("C:/Lab_projects/2016_Phylomodelling/Data/Shapefiles/states", "cb_2014_us_state_500k")
cali <- cali[cali$NAME=="California",]
cali <- spTransform(cali, proj)



# get clipped convex hull sizes for every species
cl <- makeCluster(nodes)
registerDoParallel(cl)
r <- foreach(spp = species,
             .packages=c("sp", "rgeos")) %dopar% {
                     s <- allocc[allocc$current_name_binomial==spp,]
                     h <- gConvexHull(s)
                     h <- gIntersection(h, cali)
                     if(is.null(h)) return(NA)
                     ha <- gArea(h)
                     return(ha)
             }
stopCluster(cl)

# add hull areas to existing table of species range characteristics
r <- data.frame(spp=species, hull=unlist(r))
r <- full_join(read.csv(paste0(project_stem_dir, "/git_files/data/species_occurrence_counts.csv")), r)

# ratio of hull to maxent area will determine buffer distance
r$ratio <- r$hull/r$maxent810m
hist(r$ratio)

# double square root transformation seems best for these ratios, based on some limited experimentation
trans <- function(x) sqrt(sqrt(x))
hist(trans(r$ratio))

# function to calculte buffer distance (in km)
buffer_distance <- function(hull_area, maxent_area, max_ratio, limits=c(10, 50), trans=NULL){
        fx <- function(x) x
        if(!is.null(trans)) fx <- trans
        ratio <- hull_area/maxent_area
        if(!is.finite(ratio)) ratio <- max_ratio
        if(ratio > max_ratio) ratio <- max_ratio
        dist <- fx(ratio) / fx(max_ratio) * (limits[2] - limits[1]) + limits[1]
        return(dist)
}

# check that buffer ratios are nicely distributed
buffs <- apply(r[,c("hull", "maxent810m")], 1, 
               function(x){buffer_distance(hull_area=x[1], maxent_area=x[2], 
                                           max_ratio=max(na.omit(r$ratio)), trans=trans)})
hist(buffs)


# generate and save hulls and buffers for each species
cl <- makeCluster(nodes)
registerDoParallel(cl)
md <- foreach(spp = species,
             .packages=c("sp", "rgeos")) %dopar% {
                     
                     # geospatial convex hull, clipped to coastline
                     s <- allocc[allocc$current_name_binomial==spp,]
                     h <- gConvexHull(s)
                     h <- gIntersection(h, cali)
                     
                     # areas of hull and maxent
                     ha <- gArea(h)
                     ma <- ranges$maxent810[r$spp==spp]
                     if(is.na(ma)) return(data.frame(spp=spp, hull_area=NA, hull_maxent_ratio=NA, buffer_area=NA, buffer_distance=NA))
                     
                     # buffer (calculate distance, convert to m)
                     dist <- buffer_distance(ha, ma, max(na.omit(r$ratio)), trans=trans) * 1000
                     b <- gBuffer(s, width=dist, byid=F)
                     ba <- gArea(b)
                     
                     # save hull and buffer
                     saveRDS(h, paste0(outdir, "/convex_hulls/", spp, ".rds"))
                     saveRDS(b, paste0(outdir, "/occurrence_buffers/", spp, ".rds"))
                     
                     # return some metadata (all values in meters)
                     return(data.frame(spp=spp, hull_area=ha, hull_maxent_ratio=ha/ma, buffer_area=ba, buffer_distance=dist))
             }
stopCluster(cl)


# add metatadata to master table and save
#mdd <- lapply(md, as.numeric)
mdd <- do.call("rbind", md)
r <- full_join(r, mdd)
#write.csv(r, paste0(project_stem_dir, "/git_files/data/species_occurrence_counts.csv"))



#ggplot(r, aes(hull, maxent810m, color=trans(hull_norm/mx_norm))) +
#        geom_point() +
#        geom_abline() +
#        scale_x_log10() + 
#        scale_y_log10() + 
#        scale_color_viridis()




