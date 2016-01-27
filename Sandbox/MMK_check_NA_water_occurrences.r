


# check to see how many occurrences fall in water (NA climate values)

library(raster)

# load occurrences
f <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Species/Clean_CCH", full.names=T)
f <- lapply(f, read.csv, stringsAsFactors=F)
f <- do.call("rbind", f)
f <- f[!is.na(f$longitude + f$latitude),]
coordinates(f) <- c("longitude", "latitude")
projection(f) <- '+proj=longlat +ellps=WGS84'

# load climate
r <- readRDS("C:/Lab_projects/2016_Phylomodelling/Data/Climate/BCM_normals/cwd_810m_filled3x.rds")

# extract and count
f <- spTransform(f, crs(r))
v <- extract(r, f)
length(v[is.na(v)])

# look
plot(r)
plot(f[is.na(v),], add=T)
table(f$scientificName[is.na(v)])
