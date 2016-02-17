

#### build table of species range stats, based on point, raster, and polygon data ####

# last updated February 2016 by Matthew Kling

library(raster)
library(dplyr)
library(tidyr)


####### point data #######

# load climate
r <- readRDS("C:/Lab_projects/2016_Phylomodelling/Data/Climate/BCM_normals/cwd_810m_filled3x.rds")

# load occurrences
allocc <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Species/processed3/atomic", full.names=T)
allocc <- lapply(allocc, readRDS)
allocc <- do.call("rbind", allocc)
allocc <- allocc[!is.na(allocc$longitude + allocc$latitude),]
coordinates(allocc) <- c("longitude", "latitude")
projection(allocc) <- '+proj=longlat +ellps=WGS84'
allocc <- spTransform(allocc, crs(r))

# exclude water occurrences
v <- raster::extract(r, allocc)
allocc <- allocc[!is.na(v),]

# count number of records and unique pixels per species
p <- r
values(p) <- 1:ncell(p)
p <- mask(p, r)
d <- data.frame(spp=allocc$current_name_binomial,
                px=raster::extract(p, allocc))
d <- d %>%
        group_by(spp) %>%
        mutate(nrecords=n()) %>%
        ungroup() %>%
        distinct() %>%
        group_by(spp) %>%
        summarize(ncells=n(),
                  nrecords=mean(nrecords))

# save ECDF chart of cells per species
library(ggplot2)
brks <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
p <- ggplot(d, aes(ncells)) + 
        stat_ecdf() +
        scale_x_log10(breaks=brks) +
        scale_y_continuous(breaks=seq(0,1,.1)) +
        labs(x="number of 810m pixels with occurrences",
             y="cumulative portion of species")
ggsave("C:/Lab_projects/2016_Phylomodelling/Output/Charts/pixels_per_species_ecdf.png", p, width=6, height=6, units="in")
        

# max of minimum spanning tree distances (in meters)
library(vegan)
library(rgeos)
library(sp)

getSpan <- function(species){
        z <- allocc[allocc$current_name_binomial==species,] %>%
                as.data.frame() %>%
                select(longitude, latitude) %>%
                distinct()
        coordinates(z) <- c("longitude", "latitude")
        projection(z) <- projection(allocc)
        z <- spTransform(z, crs('+proj=longlat +ellps=WGS84'))
        dst <- geosphere::distm(z)
        span <- spantree(dst)
        return(max(span$dist))
}
d$span <- sapply(d$spp, getSpan)
write.csv(d, "C:/Lab_projects/2016_Phylomodelling/git_files/data/species_max_spans.csv")
d <- select(d, -span)





####### raster and polygon data #######

model_dir <- paste0(project_stem_dir, "/Output/Maxent/v5")
richness_dir <- paste0(project_stem_dir, "/Output/Richness/V5")
pr_dir <- paste0(project_stem_dir, "/Data/Richness/point_richness_rasters")
spp_dirs <- list.dirs(model_dir, full.names=T, recursive=F)

cl <- makeCluster(nodes)
registerDoParallel(cl)
r <- foreach(spp = spp_dirs,
             .packages=c("raster", "rgeos")) %dopar% {
                     if(!file.exists(paste0(spp, "/BinaryRangePrediction.rds"))) return("no data")
                     if(!file.exists(paste0(spp, "/BufferClippedMaxent.rds"))) return("no data")
                     
                     # maxent
                     mx <- sapply(paste0(spp, "/BinaryRangePrediction", c("", "_25k", "_50k"), ".rds"), 
                                      function(x) sum(na.omit(values(readRDS(x)))))
                     
                     # buffer-clipped maxent
                     cm <- sapply(paste0(spp, "/BufferClippedMaxent", c("", "_25k", "_50k"), ".rds"), 
                                      function(x) sum(na.omit(values(readRDS(x)))))
                     
                     # convex hull
                     ch <- gArea(readRDS(paste0(project_stem_dir, "/Output/Range_polygons/convex_hulls/", basename(spp), ".rds")))
                     
                     # occurrence buffer
                     ob <- gArea(readRDS(paste0(project_stem_dir, "/Output/Range_polygons/occurrence_buffers/", basename(spp), ".rds")))
                     
                     return(c(mx, cm, ch, ob))
             }
stopCluster(cl)

good <- sapply(r, length)==8
rd <- r[good]
rd <- lapply(rd, as.vector)
rd <- do.call("rbind", rd)
rd <- cbind(basename(spp_dirs[good]), as.data.frame(rd))
names(rd) <- c("spp", "maxent810m", "maxent25km", "maxent50km", 
               "clippedMaxent810m", "clippedMaxent25km", "clippedMaxent50km", 
               "convexHull", "occurrenceBUffer")

d <- full_join(d, rd)
write.csv(d, paste0(project_stem_dir, "/git_files/data/species_occurrence_counts.csv"), row.names=F)



###### plots ######

# pairs plots, linear and log scales
png(paste0(richness_dir, "/species_range_scatterplot.png"), width=2000, height=2000)
pairs(r[,2:ncol(r)], cex=.1)
dev.off()

png(paste0(richness_dir, "/species_range_scatterplot_loglog.png"), width=2000, height=2000)
pairs(log10(r[,2:ncol(r)]), cex=.1)
dev.off()






