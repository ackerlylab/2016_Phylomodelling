
# quantify climatic niche centroid and breadth for every CA species


library(raster)
library(rgeos)
library(doParallel)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
extract <- raster::extract


cdir <- filled_climate_data_dir
spdir <- occurrence_data_dir_processed
odir <- "C:/Lab_projects/2016_Phylomodelling/Output/Niche"
species_list <- paste0(spdir, 'combined/0_Species_list_v2.rdata')
maxent_background <- paste(spdir, 'combined/10000_CA_plants_bg_810m_occbias.rdata', sep='')


# load climate data
files <- list.files(path=cdir, pattern='810m_filled3x', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt")) #use sort to ensure correct names used for clim
clim <- stack(lapply(files, function(x) readRDS(paste(cdir,x,sep="/"))))
names(clim) <-  climnames

# california boundary
cali <- rgdal::readOGR("C:/Lab_projects/2016_Phylomodelling/Data/Shapefiles/states", "cb_2014_us_state_500k")
cali <- cali[cali$NAME=="California",]
cali <- spTransform(cali, crs(readRDS(paste0(clim_dir, "/", files[1]))))


# crop climate data to california
stop("complete me")

# load occurrences
allocc <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Species/processed3/atomic", full.names=T)
allocc <- lapply(allocc, readRDS)
allocc <- do.call("rbind", allocc)
allocc <- allocc[!is.na(allocc$longitude + allocc$latitude),]
coordinates(allocc) <- c("longitude", "latitude")
projection(allocc) <- '+proj=longlat +ellps=WGS84'
allocc <- spTransform(allocc, crs(clim))

# exclude occurrences with no climate data
allocc <- allocc[!is.na(raster::extract(clim[[1]], allocc)),]
species <- unique(allocc$current_name_binomial)

# climate principal components
vals <- values(clim)
vals[,"ppt"] <- log(vals[,"ppt"]) # log-transform PPT first
vals <- scale(vals)
good <- apply(vals, 1, function(x) !is.na(sum(x)))
pca <- prcomp(vals[good,])
PC1 <- clim[[1]]
PC1[] <- NA
PC2 <- PC1
PC1[good] <- pca$x[,"PC1"]
PC2[good] <- pca$x[,"PC2"]
pc <- stack(PC1, PC2)
names(pc) <- c("pc1", "pc2")

# raster of unique pixel IDs
pixels <- pc[[1]]
pixels[] <- 1:ncell(pixels)



# compute niche stats for every species
cl <- makeCluster(nodes)
registerDoParallel(cl)
results <- foreach(spp = species,
                   .packages=c("depth", "sp", "rgeos", "raster", "dismo")) %dopar% {
                           
                           # unique occurrences of species
                           occ <- allocc[allocc$current_name_binomial==spp,]
                           occ <- occ[sample(nrow(occ), nrow(occ), replace=F),] # scramble before getting unique, to prevent systematic within-pixel directional bias
                           occ <- occ[!duplicated(extract(pixels, occ)),] # get unique
                           if(nrow(occ)<3) return("insufficient points")
                           
                           # climate values at occurrences
                           e <- extract(pc, occ)
                           
                           # centroid of occurrences (sensitive to skew)
                           occ_centroid <- apply(e, 2, mean)
                           
                           # centroid and area of convex hull of occurrences (insensitive to skew)
                           h <- e[chull(e),]
                           h <- Polygon(h)
                           occ_hull_area <- h@area
                           h <- SpatialPolygons(list(Polygons(list(h), ID="1")))
                           occ_hull_centroid <- coordinates(rgeos::gCentroid(h))
                           
                           # maxent model (corrects for sampling bias and available climate space)
                           bg <- readRDS(maxent_background)
                           mxArgs <- c("-a", "-z", "outputformat=raw", "maximumbackground=10000", "nothreshold", "nohinge")
                           mx <- maxent(pc, p=occ, a=bg, args=mxArgs)
                           pred <- predict(mx, pc)
                           
                           # center of maxent model
                           maxent_centroid <- values(pc)[which.max(pred)[1],]
                           
                           centers <- as.data.frame(rbind(occ_centroid, occ_hull_centroid, maxent_centroid))
                           centers$stat <- c("points", "hull", "maxent")
                           centers$spp <- spp
                           return(list(centers=centers, area=occ_hull_area))
                   }
stopCluster(cl)

# compile and save results
names(results) <- species
results <- results[lapply(results, class) == "list"]
saveRDS(results, paste0(odir, "/climate_niches.rds"))


##################################################
##################################################
##################################################

results <- readRDS(paste0(odir, "/climate_niches.rds"))
names(results) <- sapply(results, function(x) x$centers$spp[1])

# restrict to cali endemics
e <- read.csv("C:/Lab_projects/2016_Phylomodelling/Data/Endemism/Californian_endemic_binomials.txt", 
              stringsAsFactors=F, header=F)[,1]
results <- results[names(results) %in% e]


a <- do.call("c", lapply(results, function(x) x$area))
a <- as.data.frame(a)
a$spp <- rownames(a)

s <- do.call("rbind", lapply(results, function(x) x$centers))
s <- left_join(s, a)

g <- as.data.frame(rasterToPoints(pc))
g$stat <- "land area"


######

ggplot(rbind(select(g, pc1, stat), select(s, pc1, stat))) +
        geom_density(aes(pc1, color=stat, fill=stat), size=1, alpha=.3) +
        theme_minimal()

# land area climate space frequency heatmap
p <- ggplot(g, aes(pc1, pc2)) +
        stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
        scale_fill_viridis(option="A", trans="sqrt") +
        theme_minimal()
ggsave(paste0(odir, "/climate_space_heatmap.png"), p, width=10, height=8)

        
# make a density difference plot in env space for geog vs spp, as for geog vs sampling



###### niche breadth

ggplot(s, aes(pc1, pc2, color=sqrt(a), size=sqrt(a))) +
        geom_point() +
        scale_color_viridis(option="A") +
        facet_wrap(~stat) +
        theme_minimal()


ss <- s %>%
        mutate(pc1=plyr::round_any(pc1, 1),
               pc2=plyr::round_any(pc2, 1)) %>%
        group_by(pc1, pc2, stat) %>%
        summarize(a=mean(sqrt(a)))
p <- ggplot(ss, aes(x=pc1, y=pc2, fill=a)) +
        geom_tile() +
        scale_fill_viridis(option="A") +
        facet_wrap(~stat, ncol=2, labeller="label_both") +
        theme_minimal() +
        theme(legend.position=c(.75,.25), legend.direction="horizontal") +
        labs(title="Niche breadth as a function of niche optimum (CA endemics only)",
             fill="hull\narea    ")
ggsave(paste0(odir, "/niche_breadth_vs_optimum.png"), p, width=8, height=6)



#########################################
# CLIM vs GEOG SPACE

library(colormap)
library(grid)
library(gridExtra)


###### clim in geog space #####

colors <- colors2d(select(g, pc1, pc2))

geo <- ggplot(g, aes(x, y)) +
        geom_raster(fill=colors) +
        ggmap::theme_nothing() +
        coord_fixed()

pca <- prcomp(vals[good,])
rot <- pca$rotation
for(i in 1:ncol(rot)) rot[,i] <- rot[,i] * pca$sdev[i]^2 * 2
rot <- as.data.frame(rot)
rot$var <- rownames(rot)

env <- ggplot(g, aes(pc1, pc2)) +
        geom_point(color=colors) + 
        annotate(geom="text", x=rot$PC1, y=rot$PC2, label=rot$var, color="black", size=7) +
        annotate(geom="segment", x=0,y=0, xend=rot$PC1*.9, yend=rot$PC2*.9, color="black",
                 arrow=grid::arrow(angle=20, type="closed", length=unit(.1, "inches"))) +
        theme_minimal() +
        theme(title=element_text(size=25)) +
        labs(title="\n\n\nCalifornia climate\n\n\n")

png(paste0(odir, "/PCA_colormap.png"), width=1500, height=1000)
combo <- arrangeGrob(env, geo, ncol=2)
grid.draw(combo)
dev.off()


##### geog in clim space #######

colors <- colors2d(prcomp(select(g, x, y))$x[,1:2], colors=c("black", "red", "yellow", "cyan"))

geo <- ggplot(g, aes(x, y)) +
        geom_raster(fill=colors) +
        ggmap::theme_nothing() +
        coord_fixed()

order <- sample(nrow(g), nrow(g))
env <- ggplot(g[order,], aes(pc1, pc2)) +
        geom_point(color=colors[order]) + 
        annotate(geom="text", x=rot$PC1, y=rot$PC2, label=rot$var, color="black", size=7) +
        annotate(geom="segment", x=0,y=0, xend=rot$PC1*.9, yend=rot$PC2*.9, color="black",
                 arrow=grid::arrow(angle=20, type="closed", length=unit(.1, "inches"))) +
        theme_minimal() +
        theme(title=element_text(size=25)) +
        labs(title="\n\n\nCalifornia climate\n\n\n")

combo <- arrangeGrob(env, geo, ncol=2)
png(paste0(odir, "/PCA_colormap_2.png"), width=1500, height=1000)
grid.draw(combo)
dev.off()


##### pooled geoclim #####

library(tsne)
library(FNN)
library(MASS)
g <- select(g, -stat)
g <- scale(g)

# perform pca on all pixels
rows <- 1:nrows(g)
pc <- prcomp(g, center=T, scale=T)
cd <- pc$x[,1:2]

# perform tsne on subsample of pixels
rows <- sample(nrow(g),1000)
ts <- tsne(g[rows,], max_iter=500, perplexity=50)
cd <- prcomp(ts)$x

# nmds
rows <- sample(nrow(g),2000)
nmd <- isoMDS(dist(g[rows,]))
cd <- prcomp(nmd$points)$x

# approximate tsne scores for all pixels based on knn
nn <- knnx.index(data=as.matrix(g[rows,]), query=as.matrix(g), k=3)
#nnd <- knnx.dist(data=as.matrix(g[rows,]), query=as.matrix(g), k=3)
cd <- cd[nn[,1],]


colnames(cd) <- c("x", "y")
colors <- colors2d(cd, colors=rev(c("yellow", "red", "blue", "green")))

geo <- ggplot(as.data.frame(g), aes(x, y)) +
        geom_raster(fill=colors) +
        ggmap::theme_nothing() +
        coord_fixed(ratio=1.2)

envrows <- sample(nrow(g), nrow(g))
env <- ggplot(as.data.frame(g)[envrows,], aes(pc1, pc2)) +
        geom_point(color=colors[envrows], size=1) + 
        #annotate(geom="text", x=rot$PC1, y=rot$PC2, label=rot$var, color="black", size=7) +
        #annotate(geom="segment", x=0,y=0, xend=rot$PC1*.9, yend=rot$PC2*.9, color="black",
        #         arrow=grid::arrow(angle=20, type="closed", length=unit(.1, "inches"))) +
        theme_minimal() +
        theme(title=element_text(size=25)) +
        labs(title="\n\n\nClimate")

combo <- ggplot(as.data.frame(cd), aes(x, y)) +
        geom_point(color=colors, size=10) + 
        theme_minimal() +
        theme(title=element_text(size=25), axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) +
        labs(title="\n\nEmbedded pooled\ngeograpy & climate")


p <- arrangeGrob(arrangeGrob(combo, env, nrow=2), geo, ncol=2)
png(paste0(odir, "/PCA_colormap_tsne.png"), height=2000, width=2000)
grid.draw(p)
dev.off()









#####################################
# NICHE STATS COMPARED


# segments linking maxent and occurrence centroids in env space
s$order <- as.integer(factor(s$stat))
randspp <- s[s$spp %in% sample(unique(s$spp), 100),]

p <- ggplot() +
        geom_point(data=g, aes(pc1, pc2), color="gray75") +
        geom_path(data=randspp, aes(pc1, pc2, order=order, group=spp), color="white") +
        geom_point(data=randspp, aes(pc1, pc2, order=order, group=spp, color=stat), size=2) +
        annotate(geom="text", x=rot$PC1, y=rot$PC2, label=rot$var, color="black", size=4) +
        annotate(geom="segment", x=0,y=0, xend=rot$PC1*.8, yend=rot$PC2*.8, color="black",
                 arrow=grid::arrow(angle=20, type="closed", length=unit(.1, "inches"))) +
        theme_minimal() +
        theme(legend.position="top") +
        labs(title="Niche optima of random CA endemics:\n3 stats compared")
ggsave(paste0(odir, "/niche_definitions_scatter.png"), p, width=8, height=8)




###################
# how does species location in climate space relate to climate variable importance?
##################


# merge with variable importance data from pwds script
library(broom)
w <- readRDS(paste0(odir, "/pwds_climate_results.rds"))
w <- w[w$species %in% e,]
w$presence <- factor(w$presence)
w <- w %>%
        mutate(ppt = log10(ppt)) %>%
        mutate_each(funs(scale), cwd:ppt) %>%
        group_by(species) %>%
        do(tidy(glm(presence ~ cwd + djf + jja + ppt, data=., family="binomial"))) %>%
        filter(term != "(Intercept)") %>%
        select(species, term, estimate, p.value)
names(w)[names(w)=="species"] <- "spp"

v <- s %>%
        group_by(spp) %>%
        summarize_each(funs(mean), pc1, pc2, a) %>%
        full_join(w, .) %>%
        na.omit() %>%
        filter(abs(estimate)<5)
v$p.value[v$p.value>.2] <- .2
v$p.value[v$p.value<.0001] <- .0001

p <- ggplot(v) +
        geom_point(aes(pc1, pc2, color=estimate, size=p.value), alpha=.75) +
        scale_color_gradient2(high="darkblue", low="darkred", mid="gray80", midpoint=0) +
        scale_size(range=c(5,1), trans="log10") +
        facet_wrap(~term) +
        geom_text(data=rot, aes(x=PC1*.8, y=PC2*.8, label=var), color="black", size=3) +
        geom_segment(data=rot, aes(x=0,y=0, xend=PC1*.66, yend=PC2*.66), color="black", size=.25,
                 arrow=grid::arrow(angle=15, type="closed", length=unit(.05, "inches"))) +
        theme(panel.background=element_blank(), legend.position="top") +
        labs(title="Variable significance (from multiple logistic regression) as a function of niche optimum",
             color="Regression coefficient  ",
             size="Regression p-value  ")
ggsave(paste0(odir, "/niche_optimum_vs_variable_importance_scatter.png"), p, width=8, height=8)


