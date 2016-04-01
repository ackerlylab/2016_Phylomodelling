

library(doParallel)
library(rgeos)
library(dplyr)
library(raster)


cdir <- filled_climate_data_dir
outdir <- paste0(project_stem_dir, "/Output/Range_polygons")


# get CRS
proj <- crs(raster("C:/Lab_projects/2016_Phylomodelling/Output/Richness/V4/richness_810m_min0.tif"))

# load occurrences
allocc <- readRDS("C:/Lab_projects/2016_Phylomodelling/Data/Species/occurrences_clean.rds")
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
                     return(gArea(h))
             }
stopCluster(cl)

# get buffer distances based on hull areas
r <- data.frame(spp=species, hull=unlist(r))
r$sqrt_hull <- sqrt(r$hull) # because sqrt represents the geometric relationship between distance and area
r$buffer <- scales::rescale(r$sqrt_hull, c(10, 50))


# generate and save hulls and buffers for each species
cl <- makeCluster(nodes)
registerDoParallel(cl)
md <- foreach(spp = species,
              .packages=c("sp", "rgeos")) %dopar% {
                      
                      s <- allocc[allocc$current_name_binomial==spp,]
                      if(nrow(s)<2) return(NA)
                      dist <- r$buffer[r$spp==spp] * 1000  # species-specific buffer distance, converted to m to match projection
                      h <- gConvexHull(s) # compute convex hull
                      hb <- gBuffer(h, width=dist, byid=F)  # buffer hull before clipping, or RAM blows up
                      h <- gIntersection(h, cali)  # clip hull to CA coastline
                      pb <- gBuffer(s, width=dist, byid=F) # buffer points
                      
                      saveRDS(h, paste0(outdir, "/convex_hulls/", spp, ".rds"))
                      saveRDS(pb, paste0(outdir, "/occurrence_buffers/", spp, ".rds"))
                      saveRDS(hb, paste0(outdir, "/buffered_hulls/", spp, ".rds"))
              }
stopCluster(cl)
