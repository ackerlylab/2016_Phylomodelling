

library(raster)
library(dismo)
library(doParallel)
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)


# load 810m template
template <- raster("C:/Lab_projects/2016_Phylomodelling/Output/Richness/V5/richness_810m_min0.tif")

# load occurrences
allocc <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Species/processed3/atomic", full.names=T)
allocc <- lapply(allocc, readRDS)
allocc <- do.call("rbind", allocc)
allocc <- allocc[!is.na(allocc$longitude + allocc$latitude),]
coordinates(allocc) <- c("longitude", "latitude")
projection(allocc) <- '+proj=longlat +ellps=WGS84'
allocc <- spTransform(allocc, crs(template))

# compute records per cell
intensity <- rasterize(allocc, template, field="id", fun="count")
if(sum(values(intensity), na.rm=T) != nrow(allocc)) stop("summation error")
intensity <- mask(intensity, template)
writeRaster(intensity, "C:/Lab_projects/2016_Phylomodelling/Output/Sampling/sampling_intensity_810m.tif")
intensity <- raster("C:/Lab_projects/2016_Phylomodelling/Output/Sampling/sampling_intensity_810m.tif")

# zoom
ext <- extent(-64008.04, 101630.6, -480050.2, -337387.2)
ext <- crop(intensity, ext)

# plots
i <- as.data.frame(rasterToPoints(intensity))
e <- as.data.frame(rasterToPoints(ext))


cali_histogram_map <- function(data, varname, title, outfile, trans="linear"){
        
        names(data) <- c("x", "y", "layer")
        if(trans=="log10") data$layer <- log10(data$layer)
        
        ramp <- scale_fill_viridis()
        ramp <- scale_fill_gradientn(colours=c("black", "darkblue", "red", "yellow"))
        
        # map
        p <- ggplot(data, aes(x, y, fill=layer)) +
                geom_raster() +
                ramp +
                coord_fixed(ratio=1) +
                theme(panel.background=element_blank(), panel.grid=element_blank(),
                      axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                      legend.position="none", title=element_text(size=25)) +
                labs(title=title)
        
        # histogram data
        bw <- range(data$layer, na.rm=T)
        bw <- (max(bw) - min(bw)) / 20
        h <- data %>%
                mutate(lyr = plyr::round_any(layer, bw)) %>%
                group_by(lyr) %>%
                summarize(n=n())
        
        # histogram
        l <- ggplot(h, aes(lyr, n, fill=lyr)) + 
                geom_bar(stat="identity", width=bw) + 
                ramp +
                theme(panel.background=element_blank(), panel.grid=element_blank(),
                      axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none") +
                labs(x=varname)
        
        png(outfile, width = 1200, height = 1500)
        plot(p)
        print(l, 
              vp=viewport(x = .5, y = .6, 
                          width = unit(0.4, "npc"), height = unit(0.25, "npc"),
                          just = c("left", "bottom")))
        dev.off()
        
        
}

cali_histogram_map(i, varname="log10(records per cell)", title="CCH sampling intensity (810m)", trans="log10",
              outfile="C:/Lab_projects/2016_Phylomodelling/Output/Sampling/sampling.png")




############# environmental space ##################

# load climate data
cdir <- filled_climate_data_dir
files <- list.files(path=cdir, pattern='810m_filled3x', full.names=T)
climnames <- sort(c("cwd","djf","jja","ppt")) #use sort to ensure correct names used for predictors
s <- stack(lapply(files, readRDS))

# clip climate data to CA boundary
cali <- rgdal::readOGR("C:/Lab_projects/2016_Phylomodelling/Data/Shapefiles/states", "cb_2014_us_state_500k")
cali <- cali[cali$NAME=="California",]
cali <- spTransform(cali, crs(s))
s <- mask(s, cali)

# combine sampling and climate data
s <- stack(s, intensity)
d <- as.data.frame(rasterToPoints(s))
names(d) <- c("x", "y", climnames, "samples")
d <- d[!is.na(d$cwd),]
d$samples[is.na(d$samples)] <- 0

# transform climate variables and get first 2 PCs
d$ppt <- log10(d$ppt)
d$ppt <- scale(d$ppt)
d$cwd <- scale(d$cwd)
d$djf <- scale(d$djf)
d$jja <- scale(d$jja)
pc <- prcomp(d[,climnames])
d$pc1 <- pc$x[,1]
d$pc2 <- pc$x[,2]
ddd <- d

### three sampling methods:

# 1) area-based
d1 <- d[sample(nrow(d), 100000),]
d1$scheme <- "all cells"

# 2) exclude unsampled cells
d2 <- d[d$samples>0,]
d2 <- d2[sample(nrow(d2), 100000),]
d2$scheme <- "cells with samples"

# 3) by sampling intensity
d3 <- d[sample(nrow(d), 100000, prob=d$samples, replace=T),]
d3$scheme <- "samples"

d <- rbind(d1, d2, d3)


### pc biplot ###

png("C:/Lab_projects/2016_Phylomodelling/Output/Sampling/pc_biplot.png", width = 600, height = 600)
biplot(pc)
dev.off()


### pc colormap ###

library(colormap)
#dc <- d[d$scheme=="all cells",]
#dc <- dc[sample(nrow(dc), 10000),]
cols <- colors2d(as.matrix(ddd[,c("pc1", "pc2")]))

p <- ggplot() +
        geom_raster(data=ddd, aes(x, y), fill=cols) +
        coord_fixed(ratio=1) +
        theme(panel.background=element_blank(), panel.grid=element_blank(),
              axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
              legend.position="none", title=element_text(size=25)) +
        labs(title="California climate")

l <- ggplot() +
        geom_point(data=ddd, aes(pc1, pc2), color=cols) +
        theme(panel.background=element_blank(), panel.grid=element_blank(),
              #axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
              legend.position="none", title=element_text(size=25))

png("C:/Lab_projects/2016_Phylomodelling/Output/Sampling/pc_map.png", width = 1200, height = 1500)
plot(p)
print(l, 
      vp=viewport(x = .5, y = .6, 
                  width = unit(0.4, "npc"), height = unit(0.3, "npc"),
                  just = c("left", "bottom")))
dev.off()


### heatmap ###

p <- ggplot(d, aes(pc1, pc2)) +
        stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
        scale_fill_gradientn(colours=c("black", "blue", "red", "yellow")) +
        facet_wrap(~scheme) +
        theme_minimal() +
        theme(legend.position="top", panel.grid=element_blank()) +
        labs(fill="frequency", title="CCH sampling effort in environmental space")
ggsave("C:/Lab_projects/2016_Phylomodelling/Output/Sampling/env_heatmap.png", p, width=12, height=6, units="in")



### difference heatmap ###

library(MASS)
lims <- c(range(d$pc1), range(d$pc2))
k1 <- kde2d(d$pc1[d$scheme=="all cells"], d$pc2[d$scheme=="all cells"], lims=lims, n=100)
k2 <- kde2d(d$pc1[d$scheme=="cells with samples"], d$pc2[d$scheme=="cells with samples"], lims=lims, n=100)
k3 <- kde2d(d$pc1[d$scheme=="samples"], d$pc2[d$scheme=="samples"], lims=lims, n=100)

flatten <- function(f){
        df <- as.data.frame(f$z)
        names(df) <- paste0("y", f$y)
        df$x <- f$x
        df <- tidyr::gather(df, y, z, -x)
        df$y <- as.numeric(sub("y", "", df$y))
        df
}

k1 <- flatten(k1)
k2 <- flatten(k2)
k3 <- flatten(k3)

k3$diff <- (k3$z - k1$z)
p <- ggplot(k1, aes(x, y, fill=diff)) + 
        geom_raster() + 
        scale_fill_gradient2(high="green", mid="black", low="red", midpoint=0) +
        #scale_fill_gradient(high="yellow", low="black") +
        theme_minimal() +
        theme(legend.position="top", panel.grid=element_blank()) +
        labs(fill="sample density minus true density", 
             title="CCH sampling bias in environmental space",
             x="pc1", y="pc2")
ggsave("C:/Lab_projects/2016_Phylomodelling/Output/Sampling/env_heatmap_bias.png", p, width=6, height=7, units="in")



### contours ###

contourLevel <- function(x,y,prob=0.95) {
        kk <- MASS::kde2d(x,y)
        dx <- diff(kk$x[1:2])
        dy <- diff(kk$y[1:2])
        sz <- sort(kk$z)
        c1 <- cumsum(sz) * dx * dy
        approx(c1, sz, xout = 1 - prob)$y
}

c1 <- contourLevel(d1$pc1, d1$pc2, prob=c(.5, .9, .99))
c2 <- contourLevel(d2$pc1, d2$pc2, prob=c(.5, .9, .99))
c3 <- contourLevel(d3$pc1, d3$pc2, prob=c(.5, .9, .99))


brk1 <- c(c1[1], c2[1], c3[1])
contours1 <- mapply(function(data, b) stat_density2d(data=data, 
                                                     aes(x=pc1, y=pc2, color=scheme, fill=scheme), 
                                                     breaks=b, geom="polygon", alpha=.1, size=.1), 
                    plyr::dlply(d, plyr::.(scheme)), brk1)

brk2 <- c(c1[2], c2[2], c3[2])
contours2 <- mapply(function(data, b) stat_density2d(data=data, 
                                                     aes(x=pc1, y=pc2, color=scheme, fill=scheme), 
                                                     breaks=b, geom="polygon", alpha=.1, size=.1), 
                    plyr::dlply(d, plyr::.(scheme)), brk2)

brk3 <- c(c1[3], c2[3], c3[3])
contours3 <- mapply(function(data, b) stat_density2d(data=data, 
                                                     aes(x=pc1, y=pc2, color=scheme, fill=scheme), 
                                                     breaks=b, geom="polygon", alpha=.1, size=.1), 
                    plyr::dlply(d, plyr::.(scheme)), brk3)


p <- ggplot() + 
        contours3 + contours2 + contours1 +
        theme_minimal() +
        scale_x_continuous(expand=c(0.1,0)) +
        scale_y_continuous(expand=c(0.1,0)) +
        facet_wrap(~scheme) +
        labs(title="Sampling intensity in environmental space:\n(50%, 90%, and 99% contours for 3 distributions)") +
        theme(legend.position="top")
ggsave("C:/Lab_projects/2016_Phylomodelling/Output/Sampling/env_contours.png", p, width=12, height=6, units="in")




#########

brk <- c()

# compute probability contour level for each month
for(month in unique(d$month)){
        dd <- d[d$variable==variable & d$month==month,]
        dd <- dd[is.finite(dd$intercept) & is.finite(dd$trendp),]
        brk <- c(brk, contourLevel(dd$intercept, dd$trendp, prob=prob))
}

# generate graph series with month-specific breaks
#library(plyr)
contours <- mapply(function(data, b) stat_density2d(data=data, 
                                                    aes(x=intercept, y=trendp, color=month, fill=month), 
                                                    breaks=b, geom="polygon", alpha=.2, size=.1), 
                   dlply(dv, .(month)), brk)

p <- ggplot(d) +
        geom_hline(yintercept=0, color="gray") +
        contours +
        geom_smooth(data=dv, aes(x=intercept, y=trendp, color=month), method=lm, se=F) +
        geom_polygon(data=pv, aes(mean, shift, order=as.integer(month)), fill=NA, color="black", linetype=2) +
        geom_smooth(data=pv, aes(mean, shift), method=lm, se=F, color="black") +
        geom_point(data=pv, aes(mean, shift, color=month), size=3) +
        labs(x="intercept", y="trend", title=paste("Seasonal and spatial dynamics of recent change in", translate(variable, "words"))) +
        whiteness()
