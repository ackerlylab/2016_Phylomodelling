

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

# get distribution of clipped convex hull areas
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

# get buffer distances based on hull areas
r <- data.frame(spp=species, hull=unlist(r))
r$sqrt_hull <- sqrt(r$hull)
r$buffer <- scales::rescale(r$sqrt_hull, c(10, 50))


# generate and save hulls and buffers for each species
cl <- makeCluster(nodes)
registerDoParallel(cl)
md <- foreach(spp = species,
              .packages=c("sp", "rgeos")) %dopar% {
                      
                      # geospatial convex hull, clipped to coastline
                      s <- allocc[allocc$current_name_binomial==spp,]
                      h <- gConvexHull(s)
                      h <- gIntersection(h, cali)
                      if(is.null(h)) return(NA)
                      
                      # hull area
                      ha <- gArea(h)
                      
                      # buffer (calculate distance, convert to m)
                      dist <- r$buffer[r$spp==spp]
                      b <- gBuffer(s, width=dist, byid=F)
                      
                      # save hull and buffer
                      saveRDS(h, paste0(outdir, "/convex_hulls/", spp, ".rds"))
                      saveRDS(b, paste0(outdir, "/occurrence_buffers/", spp, ".rds"))
              }
stopCluster(cl)
