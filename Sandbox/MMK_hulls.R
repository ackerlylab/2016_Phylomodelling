

library(raster)
library(doParallel)
library(rgeos)
library(dplyr)

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

# clip occurrences to california boundary
cali <- rgdal::readOGR("C:/Lab_projects/2016_Phylomodelling/Data/Shapefiles/states", "cb_2014_us_state_500k")
cali <- cali[cali$NAME=="California",]
cali <- spTransform(cali, proj)
ao <- over(allocc, cali)
allocc <- allocc[!is.na(ao$GEOID),]

# load table of maxent range sizes
ranges <- read.csv(paste0(project_stem_dir, "/git_files/data/species_occurrence_counts.csv"))

species <- unique(allocc$current_name_binomial)


cl <- makeCluster(7)
registerDoParallel(cl)
r <- foreach(spp = species,
             .packages=c("sp", "rgeos")) %dopar% {
                     
                     # compute geospatial convex hull and clip to coastline
                     s <- allocc[allocc$current_name_binomial==spp,]
                     h <- gConvexHull(s)
                     h <- gIntersection(h, cali)
                     
                     ha <- gArea(h)
                     #ma <- ranges$maxent810[ranges$spp==spp]
                     #distance <- ha / ma
                     
                     return(ha)
             }
stopCluster(cl)

r <- rr$hull
r <- data.frame(spp=species, hull=r)
r <- full_join(r, ranges)
r$hull_norm <- scales::rescale(r$hull)
r$mx_norm <- scales::rescale(r$maxent810m)

ggplot(r, aes(hull_norm, mx_norm, color=log10(hull_norm/mx_norm))) +
        geom_point() +
        geom_abline() +
        scale_x_log10() + 
        scale_y_log10() + 
        scale_color_viridis()

r$ratio <- sqrt(sqrt(r$hull/r$maxent810m))
hist(r$ratio)

rng <- range(r$ratio[is.finite(r$ratio)])


buffer_distance <- function(hull_area, maxent_area, limits=c(10, 50)){
        sqrt(sqrt(log10(ha/ma))) / rng[2] * (buff_lims[2] - buff_lims[1]) + buff_lims[1]
}


# buffer limits, in km


cl <- makeCluster(7)
registerDoParallel(cl)
r <- foreach(spp = species,
             .packages=c("sp", "rgeos")) %dopar% {
                     
                     # compute geospatial convex hull and clip to coastline
                     s <- allocc[allocc$current_name_binomial==spp,]
                     h <- gConvexHull(s)
                     h <- gIntersection(h, cali)
                     
                     ha <- gArea(h)
                     ma <- ranges$maxent810[ranges$spp==spp]
                     
                     # compute buffer distance based on hull/maxent ratio
                     dist <- 
                     dist <- dist
                     
                     bp <- gBuffer(s, dist)
                     return(ha)
             }
stopCluster(cl)
