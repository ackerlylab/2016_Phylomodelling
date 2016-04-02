



library(raster)
library(rgeos)
library(doParallel)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(colormap)
library(grid)
library(gridExtra)

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
cali <- spTransform(cali, crs(readRDS(paste0(cdir, "/", files[1]))))

clim <- mask(clim, cali)



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




g <- as.data.frame(rasterToPoints(pc))
g$stat <- "land area"



#########################################
# CLIM vs GEOG SPACE


###### clim in geog space #####

colors <- colors2d(select(g, pc1, pc2), c("yellow", "red", "black", "cyan"))

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
        theme(plot.title=element_text(size=35)) +
        labs(title="\n\n\nCalifornia climate\n\n\n")

png(paste0(odir, "/PCA_colormap.png"), width=1500, height=1000)
combo <- arrangeGrob(env, geo, ncol=2)
grid.draw(combo)
dev.off()


##### geog in clim space #######

colors <- colors2d(prcomp(select(g, x, y))$x[,1:2], colors=c("yellow", "red", "black", "cyan"))

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
        theme(plot.title=element_text(size=35)) +
        labs(title="\n\n\nCalifornia climate\n\n\n")

combo <- arrangeGrob(env, geo, ncol=2)
png(paste0(odir, "/PCA_colormap_2.png"), width=1500, height=1000)
grid.draw(combo)
dev.off()


##### pooled geoclim #####

library(tsne)
library(FNN)
library(MASS)
g <- dplyr::select(g, -stat)
g <- scale(g) # this may inappropriately give equal weight to longitude and latitude


for(method in c("pca", "tsne", "nmds", "som", "sammon")){
        
        if(method=="pca"){
                pc <- prcomp(g, center=T, scale=T)
                cd <- pc$x[,1:2]
                pointsize=.5
        }
        
        if(method=="nmds"){
                rows <- sample(nrow(g),1000)
                nmd <- isoMDS(dist(g[rows,]))
                cd <- prcomp(nmd$points)$x
                pointsize=10
        }
        
        if(method=="sammon"){
                rows <- sample(nrow(g),3000)
                sam <- sammon(dist(g[rows,]))
                cd <- sam$points
                pointsize=7
        }
        
        if(method=="som"){
                require(kohonen)
                n <- 25
                rows <- sample(nrow(g),10000)
                somdata <- g[rows,]
                fit <- somgrid(n,n)
                fit <- som(somdata, grid=fit, rlen=500)
                cd <- fit$grid$pts[fit$unit.classif,]
        }
        
        if(method=="tsne"){
                n <- 1000
                rows <- sample(nrow(g), n)
                ts <- tsne(g[rows,], max_iter=500, perplexity=n/10, whiten=F)
                cd <- prcomp(ts)$x
                pointsize=10
        }
        
        
        # approximate scores for all pixels based on knn
        if(method != "pca"){
                nn <- knnx.index(data=as.matrix(g[rows,]), query=as.matrix(g), k=1)
                cd <- cd[nn[,1],]
        }
        
        cd <- as.data.frame(cd)
        colnames(cd) <- c("x", "y")
        
        if(method=="som") cd <- cd %>% group_by(x, y) %>% mutate(n=n())
        
        
        #colors <- colors2d(cd[,c("x", "y")], colors=c("yellow", "red", "black", "cyan"))
        colors <- colors2d(cd[,c("x", "y")], colors=c("red", "yellow", "cyan", "black"))
        
        style <- theme(plot.title=element_text(size=45), axis.title=element_blank(), 
                       axis.ticks=element_blank(), axis.text=element_blank(),
                       panel.grid=element_blank(), panel.background=element_blank())
        
        geo <- ggplot(as.data.frame(g), aes(x, y)) +
                geom_raster(fill=colors) +
                coord_fixed(ratio=1.2) +
                labs(title="\nGeographic space") +
                style
        
        envrows <- sample(nrow(g), nrow(g)) # scramble order to prevent overplotting bias
        env <- ggplot(as.data.frame(g)[envrows,], aes(-pc2, -pc1)) +
                geom_point(color=colors[envrows], size=1) + 
                labs(title="\nClimate space") +
                coord_fixed() +
                style
        
        if(method=="som"){
                combo <- ggplot(cd, aes(x, -y)) +
                        geom_raster(fill=colors, alpha=scales::rescale(ecdf(cd$n)(cd$n))) + 
                        labs(title=paste0("\nEco space?\n(pooled ", toupper(method), ")")) +
                        coord_fixed() +
                        style
        } else{
                combo <- ggplot(cd, aes(y, -x)) +
                        geom_point(color=colors, size=pointsize) + 
                        labs(title=paste0("\nEco space?\n(pooled ", toupper(method), ")")) +
                        #coord_fixed() +
                        style
        }
        
        p <- arrangeGrob(env, combo, geo, ncol=3)
        png(paste0(odir, "/ecospace_", method, ".png"), width=3000, height=1000)
        grid.draw(p)
        dev.off()
        
}



##############################
# how regionally and locally rare is each pixels climate?

library(MASS)
library(FNN)


#### regional #####

# kernel density in climate space
vals <- na.omit(values(pc))
vals <- apply(vals, 2, scales::rescale)
kdd <- kde2d(vals[,1], vals[,2], n=100, h=.5)
kd <- kdd$z
#image(kd); points(vals[sample(nrow(vals), 5000),], cex=.2)

# kd row and col numbers for each value
dens <- vals
dens[,1] <- knnx.index(data=kdd$x, query=dens[,1], k=1)
dens[,2] <- knnx.index(data=kdd$y, query=dens[,2], k=1)
dens <- apply(dens, 1, function(x) kd[x[1],x[2]])

density <- pc[[1]]
density[!is.na(values(density))] <- dens
names(density) <- "density"

gd <- as.data.frame(rasterToPoints(density))
geo <- ggplot(gd, aes(x, y, fill=-density+max(gd$density))) + 
        geom_raster() + 
        scale_fill_viridis() +
        labs(title="regionwide rarity of local climate\n",
             fill="rarity") +
        coord_fixed() +
        style +
        theme(legend.position="top") +
        guides(fill=guide_colourbar(barwidth=30))

ed <- cbind(as.data.frame(vals), dens)
env <- ggplot(ed, aes(pc1, pc2, color=-dens+max(gd$density))) +
        geom_point() + 
        scale_color_viridis() +
        style +
        theme(rect=element_blank(), legend.position="none")

png(paste0(odir, "/climate_rarity.png"), width=1200, height=1500)
plot(geo)
print(env, 
      vp=viewport(x = .45, y = .55, 
                  width = unit(0.5, "npc"), height = unit(0.35, "npc"),
                  just = c("left", "bottom")))
dev.off()






################ local ###############


focal_stack <- function(r, fun, radius, circular=T, ...){
        
        # todo:
        # allow for arbitrary number of output layers
        # implement parallel processing
        
        rad <- radius
        r <- extend(r, rad)
        a <- as.array(r)
        atemplate <- a
        rtemplate <- r
        
        #nchunks <- 6
        #groups <- as.integer(cut(1:nrow(a), seq(0, nrow(a), nrow(a)/nchunks)))
        
        # circular mask
        if(circular){
                diam <- rad * 2 + 1
                m <- matrix(1, diam,diam)
                m <- cbind(rep(1:nrow(m), ncol(m)), rep(1:ncol(m), each=nrow(m)))
                m <- m - m[nrow(m)/2+.5,]
                m <- apply(m, 1, function(x) sqrt(sum(x^2)))
                m <- matrix(m, diam, diam)
                m[m <= m[rad+1, 1]] <- 1
                m[m != 1] <- NA
                m <- array(rep(m, dim(a)[3]), c(diam, diam, dim(a)[3]))
        }
        
        for(i in 1:nrow(a)){
                for(j in 1:ncol(a)){
                        if(is.na(a[i,j,1])) next()
                        w <- a[(i-rad):(i+rad), (j-rad):(j+rad),]
                        if(circular) w <- w * m
                        atemplate[i,j,] <- fun(w, ...)
                }
        }
        
        values(rtemplate) <- atemplate
        return(rtemplate)
}

isolation <- function(w){
        require(depth)
        center <- nrow(w)/2+.5
        e <- w[center, center,]
        w <- apply(w, 3, function(d)d)
        w <- na.omit(w)
        i <- .5 - depth(e, w)
        return(i)
}

heterogeneity <- function(w){
        center <- nrow(w)/2+.5
        w <- apply(w, 3, function(d)d)
        w <- na.omit(w)
        h <- chull(w)
        h <- c(h, h[1])
        h <- w[h,]
        h <- Polygon(h, hole=F)
        h <- h@area
        return(h)
}


radius <- 100 / .810 # 100km search radius, in pixels

iso <- focal_stack(pc, isolation, radius)[[1]]
writeRaster(iso, "C:/Lab_projects/2016_Phylomodelling/Output/Climate/isolation.tif", overwrite=T)

het <- focal_stack(pc, heterogeneity, radius)[[1]]
writeRaster(het, "C:/Lab_projects/2016_Phylomodelling/Output/Climate/heterogeneity.tif", overwrite=T)

locrar <- iso * het
writeRaster(locrar, "C:/Lab_projects/2016_Phylomodelling/Output/Climate/local_rarity.tif", overwrite=T)




