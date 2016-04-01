

# prepare a table of clean final occurrence points

library(rgeos)
library(dplyr)
library(raster)

cdir <- filled_climate_data_dir
outdir <- paste0(project_stem_dir, "/Output/Range_polygons")


# get CRS
proj <- crs(raster("C:/Lab_projects/2016_Phylomodelling/Output/Richness/V4/richness_810m_min0.tif"))

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

# save
saveRDS(allocc, "C:/Lab_projects/2016_Phylomodelling/Data/Species/occurrences_clean.rds")