
# generate raster representing distance to nearest occurrence

# last updated February 2016 by Matthew Kling

library(raster)
library(dplyr)
library(tidyr)
library(doParallel)
library(ggplot2)
extract <- raster::extract

######### file paths ##########
# source user_parameters.r before running #

clim_dir <- filled_climate_data_dir
model_dir <- paste0(project_stem_dir, "/Output/Maxent/v5")
richness_dir <- paste0(project_stem_dir, "/Output/Richness/V5")
pr_dir <- paste0(project_stem_dir, "/Data/Richness/point_richness_rasters")


######### climate data ###########
files <- list.files(path=clim_dir, pattern='810m_filled3x', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt")) #use sort to ensure correct names used for predictors


######### maxent models ###############
spp_dirs <- list.dirs(model_dir, full.names=T, recursive=F)


####### point data #######

# load climate raster for reference
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

stop()

############## experiment with a single species ###############

spp <- spp_dirs[grepl("Adenophyllum cooperi", spp_dirs)]

# maxent prediction surface
m <- readRDS(paste0(spp, "/ModelObject.rdata"))
p <- stack(lapply(files,function(x) readRDS(paste(clim_dir,x,sep="/"))))
names(p) <- climnames
pred <- predict(m, p)

# apply multiplier based on gaussian kernel distance from occurrence points
occ <- allocc[allocc$current_name_binomial==basename(spp),]
dst <- mask(distanceFromPoints(r, occ), r)
#gauss <- function(x, sigma) 1 / (sqrt(2*pi)) * exp(-.5 * (x/sigma) ^ 2)
gauss <- function(x, sigma) exp(-.5 * (x/sigma) ^ 2)
sigma <- 500 # distance in km
mult <- calc(dst, function(x) gauss(x, sigma*1000))
combo <- pred * mult


set.seed(123)
bg <- allocc[sample(nrow(allocc), 1000),]

getBinary <- function(x){
        eval <- evaluate(extract(x, occ), extract(x, bg))
        thresh <- threshold(eval, stat="spec_sens")
        x <- reclassify(x, c(-Inf, thresh, 0, thresh, Inf, 1))
        return(x)
}

bpred <- getBinary(pred)
bcombo <- getBinary(combo)

# plot
s <- as.data.frame(rasterToPoints(stack(pred, mult, combo, bpred, bcombo)))
s <- s[sample(nrow(s), 1000),]
names(s) <- c("x", "y", "pred", "mult", "combo", "bpred", "bcombo")
ggplot(s, aes(pred, mult, color=paste(bpred, bcombo))) + 
        geom_point(size=6) + 
        #viridis::scale_color_viridis() + 
        theme_minimal() +
        labs(x="maxent", y="gaussian distance")








################## deprecated attempt1 is below this line #####################

############# what's the right kernel distance? get stats from many species ##############

# background points
set.seed(123)
bg <- allocc[sample(nrow(allocc), 1000),]

cl <- makeCluster(nodes)
registerDoParallel(cl)
results <- foreach(spp = spp_dirs[sample(length(spp_dirs), 200)],
                   .packages=c("raster", "dismo")) %dopar% {
                           
                           # function to calculate and apply threshold
                           getBinary <- function(x){
                                   eval <- evaluate(extract(x, occ), extract(x, bg))
                                   thresh <- threshold(eval, stat="spec_sens")
                                   x <- reclassify(x, c(-Inf, thresh, 0, thresh, Inf, 1))
                                   return(x)
                           }
                           
                           # occurrences
                           occ <- allocc[allocc$current_name_binomial==basename(spp),]
                           
                           # maxent prediction
                           m <- readRDS(paste0(spp, "/ModelObject.rdata"))
                           p <- stack(lapply(files,function(x) readRDS(paste(clim_dir,x,sep="/"))))
                           names(p) <- climnames
                           pred <- predict(m, p)
                           bpred <- getBinary(pred)
                           
                           # apply multiplier based on gaussian kernel distance from occurrence points
                           dst <- mask(distanceFromPoints(r, occ), r)
                           #gauss <- function(x, sigma) 1 / (sqrt(2*pi)) * exp(-.5 * (x/sigma) ^ 2)
                           gauss <- function(x, sigma) exp(-.5 * (x/sigma) ^ 2)
                           
                           # test various distances
                           sigmas <- c(10, 20, 50, 100, 200, 500)
                           f <- data.frame(Var1=NULL, Freq=NULL, sigma=NULL)
                           for(sigma in sigmas){
                                   mult <- calc(dst, function(x) gauss(x, sigma*1000))
                                   combo <- pred * mult
                                   bcombo <- getBinary(combo)
                                   fs <- as.data.frame(table(na.omit(values(bpred + bcombo*10))))
                                   fs$sigma <- sigma
                                   f <- rbind(f, fs)
                           }
                           
                           f$spp <- basename(spp)
                           return(f)
                   }
stopCluster(cl)


d <- do.call("rbind", results)
all <- expand.grid(Var1=unique(d$Var1), sigma=unique(d$sigma), spp=unique(d$spp))
d <- full_join(d, all)
d$Freq[is.na(d$Freq)] <- 0
sigmas <- c(10, 20, 50, 100, 200, 500)

ggplot(d[d$spp %in% weird,]) +
        geom_line(aes(sigma, Freq, group=spp)) +
        #geom_smooth(aes(sigma, Freq)) +
        facet_wrap(~Var1, scales="free") +
        scale_x_log10(breaks=sigmas)


s <- d %>%
        mutate(Var1=as.integer(as.character(Var1))) %>%
        filter(Var1==10) %>%
        group_by(spp) %>%
        summarize(sigmax=sigma[Freq==max(Freq)][1])

weird <- s$spp[s$sigmax==500]


# goals:
# minimize climatically unsuitable land that is added simply due to proximity (category 10)


# want to minimize the proportion 
s <- d %>%
        mutate(Var1=paste0("p", Var1)) %>%
        spread(Var1, Freq) %>%
        mutate(prop = p10 / (p1 + p11)) #### should this be p11, or p1 + p11??

ggplot(s, aes(p1, fill=factor(sigma), color=factor(sigma))) + 
        geom_density(alpha=.2) +
        scale_x_log10()

s <- s %>%
        group_by(sigma) %>%
        summarize(prop=mean(prop))


