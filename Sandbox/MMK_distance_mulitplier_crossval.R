



# set directories
cdir <- filled_climate_data_dir
spdir <- occurrence_data_dir_processed
odir <- maxent_output_dir
species_list <- paste0(spdir, 'combined/0_Species_list_v2.rdata')
maxent_background <- paste(spdir, 'combined/10000_CA_plants_bg_810m_occbias.rdata', sep='')

# load libraries
options(java.parameters = "-Xmx1g" )
require(dismo)
require(rJava)
require(raster) 

# load climate data
files <- list.files(path=cdir, pattern='810m_filled3x', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt")) #use sort to ensure correct names used for predictors
predictors <- stack(lapply(files,function(x) readRDS(paste(cdir,x,sep="/"))))
names(predictors) <-  climnames 

# load occurrences
allocc <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Species/processed3/atomic", full.names=T)
allocc <- lapply(allocc, readRDS)
allocc <- do.call("rbind", allocc)
allocc <- allocc[!is.na(allocc$longitude + allocc$latitude),]
coordinates(allocc) <- c("longitude", "latitude")
projection(allocc) <- '+proj=longlat +ellps=WGS84'
allocc <- spTransform(allocc, crs(predictors))

# exclude occurrences with no climate data
allocc <- allocc[!is.na(raster::extract(predictors[[1]], allocc)),]
species <- unique(allocc$current_name_binomial)

# raster of unique pixel IDs
pixels <- predictors[[1]]
pixels[] <- 1:ncell(pixels)


# function to split points into four quadrants
### but wait -- while normally we want spatial independence to fight autocorrelation, that's not true with the distance method.
makeFolds <- function(points){
        crds <- coordinates(points)
        lat <- crds[,1]>median(crds[,1])
        lonN <- !lat & crds[,2]>median(crds[,2][!lat])
        lonS <- lat & crds[,2]>median(crds[,2][lat])
        f <- as.integer(lat)
        f[!lat] <- as.integer(lonN[!lat])
        f[lat] <- as.integer(lonS[lat]) + 2
        return(f)
}


spp <- species[1000]

# unique occurrences of species
occ <- allocc[allocc$current_name_binomial==spp,]


dst <- mask(distanceFromPoints(predictors[[1]], occ), predictors[[1]])
#gauss <- function(x, sigma) 1 / (sqrt(2*pi)) * exp(-.5 * (x/sigma) ^ 2)
gauss <- function(x, sigma) exp(-.5 * (x/sigma) ^ 2)
sigma <- 500 # distance in km
mult <- calc(dst, function(x) gauss(x, sigma*1000))


px <- extract(pixels, occur)
occur <- occur[!duplicated(px),]

# some maxent prep
mxArgs <- c("-a", "-z", "outputformat=raw", "maximumbackground=10000", "nothreshold", "nohinge")
bg <- readRDS(maxent_background)

occur$fold <- makeFolds(occur)
for(fold in unique(occur$fold)){
        train <- occur[occur$fold != fold,]
        test <- occur[occur$fold == fold,]
        
        # fit maxent model
        mx <- try(maxent(predictors, p=train, a=bg, args=mxArgs))
        
        
}



#mx.dir <- paste(odir, maxent_run_version, spp, sep="/")
#if(file.exists(mx.dir) == F) {dir.create(mx.dir, recursive=T)}

#saveRDS(mx, file = paste(mx.dir, 'ModelObject.rdata', sep='/'))




##### cross-eval to test proper buffer distance

for(spp in spp_dirs){
        for(fold in 1:4){
                
                # separate training and evaluation points
                
                
        }
}
