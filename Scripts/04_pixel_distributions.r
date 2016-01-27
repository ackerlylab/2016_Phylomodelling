


# check to see how many occurrences fall in water (NA climate values)
# and distribution of number of unique grid cells with records per species

library(raster)
library(dplyr)
library(tidyr)

# load climate
r <- readRDS("C:/Lab_projects/2016_Phylomodelling/Data/Climate/BCM_normals/cwd_810m_filled3x.rds")

# load occurrences
f <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Species/processed3", full.names=T)
f <- lapply(f[3:length(f)], readRDS) # the first two files aren't species records
f <- do.call("rbind", f)
f <- f[!is.na(f$longitude + f$latitude),]
coordinates(f) <- c("longitude", "latitude")
projection(f) <- '+proj=longlat +ellps=WGS84'
f <- spTransform(f, crs(r))


######## WATER OCCURRENCES #######

# extract and count
v <- raster::extract(r, f)
#length(v[is.na(v)])

# look
#plot(r)
#plot(f[is.na(v),], add=T)
#table(f$current_name_binomial[is.na(v)])


####### PIXELS/SPECIES FREQUENCY #######

f <- f[!is.na(v),]
p <- r
values(p) <- 1:ncell(p)
p <- mask(p, r)
d <- data.frame(spp=f$current_name_binomial,
                px=raster::extract(p, f))
d <- d %>%
        group_by(spp) %>%
        mutate(nrecords=n()) %>%
        ungroup() %>%
        distinct() %>%
        group_by(spp) %>%
        summarize(ncells=n(),
                  nrecords=mean(nrecords))

write.csv(d, "C:/Lab_projects/2016_Phylomodelling/git_files/data/species_occurrence_counts.csv")

library(ggplot2)
brks <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
p <- ggplot(d, aes(ncells)) + 
        stat_ecdf() +
        scale_x_log10(breaks=brks) +
        scale_y_continuous(breaks=seq(0,1,.1)) +
        labs(x="number of 810m pixels with occurrences",
             y="cumulative portion of species")
ggsave("C:/Lab_projects/2016_Phylomodelling/Output/Charts/pixels_per_species_ecdf.png", p, width=6, height=6, units="in")
        
dd <- d %>%
        group_by(ncells) %>%
        summarize(nspp=n()) %>%
        arrange(ncells) %>%
        as.data.frame()

p <- ggplot(dd, aes(ncells, nspp)) +
        geom_point() +
        scale_x_log10(breaks=brks) + 
        scale_y_log10(breaks=brks)




