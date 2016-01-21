

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(scales)
library(ggmap)

# input data
spp_dirs <- list.dirs("E:/Phylo_modelling/Output/Maxent/V2", recursive=F)
cdir <- 'E:/BCM/CA_2014/Summary/HST/Normals_30years/'


# compile a master data frame of all points
openData <- function(x){
        file <- paste0(x, "/Suitability_data.rdata")
        if(file.exists(file)){
                x <- readRDS(paste0(x, "/Suitability_data.rdata"))
                as.data.frame(x$occur)
        }
}
d <- do.call("rbind", lapply(spp_dirs, openData))
names(d)[names(d)=="current_name_binomial"] <- "spp"


# quick peek at some distirbutions
hist(d$suitability)
hist(log(d$prop_max.suit.occ))
hist(d$prop_max.suit.all)
#plot(d$prop_max.suit.all, d$prop_max.suit.occ)


# for various thresholds, summarize outlier frequency by species
d <- d %>%
        group_by(spp) %>%
        mutate(prop_max = suitability / max(suitability, na.rm=T))
s <- d %>%
        filter(!is.na(suitability)) %>%
        #select(-id) %>%
        distinct(longitude, latitude, cwd, djf, jja, ppt) %>%
        summarize(records = n(),
                  t0001 = length(prop_max[prop_max<=.001]) / n(),
                  t0003 = length(prop_max[prop_max<=.003]) / n(),
                  t0010 = length(prop_max[prop_max<=.01]) / n(),
                  t0030 = length(prop_max[prop_max<=.03]) / n(),
                  t0100 = length(prop_max[prop_max<=.1]) / n(),
                  t0300 = length(prop_max[prop_max<=.3]) / n(),
                  t1000 = length(prop_max[prop_max<=1]) / n()) %>%
        gather(threshold, anomaly_rate, -spp, -records) %>%
        mutate(threshold = as.numeric(paste0(substr(threshold, 2,2), ".", substr(threshold, 3,5))))

p <- ggplot(d, aes(prop_max)) + 
        stat_ecdf(size=1) + 
        scale_x_log10(breaks=10^(-10:1), label=scientific_format(digits=1)) + 
        scale_y_log10(breaks=10^(-10:1), label=scientific_format(digits=1)) +
        labs(x="threshold (prop max obs suitability)", 
             y="proportion of total records flagged",
             title="Distribution of maxent outliers, overall")
ggsave("E:/Phylo_modelling/Output/Charts/outlier_curve_overall_V2.png", p, width=8, height=6)

p <- ggplot(s, aes(anomaly_rate, color=factor(threshold))) + 
        stat_ecdf(size=1) +
        labs(x="proportion of within-species occurrences <= threshold",
             y="cumulative proportion of species",
             color="threshold\n(prop.max.occ)",
             title="Distirbution of maxent outliers, species-wise") +
        scale_x_log10(breaks=10^(-3:1))
ggsave("E:/Phylo_modelling/Output/Charts/outlier_curves_V2.png", p, width=8, height=6)




# PCA of all CA climate values
cdir <- 'E:/BCM/CA_2014/Summary/HST/Normals_30years/'
cfiles <- list.files(cdir, pattern="filled", full.names=T)
clim <- lapply(cfiles, readRDS)
clim <- do.call("stack", clim)
clim <- na.omit(values(clim))
colnames(clim) <- substr(basename(cfiles), 1, 3)
pc <- prcomp(clim, center=T, scale=T)


# define some spatial stuff for later use
prj <- projection(readRDS(paste0(spp_dirs[1], "/Suitability_data.rdata"))$bg)
newprj <- '+proj=longlat +ellps=WGS84'
map <- ggmap(get_map("california", source="stamen", maptype="terrain-background", zoom=6))
grid <- readRDS(paste0(cdir, "/background50by50km.Rdata"))


# point maps of maxent outliers in geographic and climate space
outlier_map <- function(x){
        
        if(!is.null(dev.list())) tryyyyyyy <- try(dev.off())
        
        require(ggplot2)
        require(gridExtra)
        require(raster)
        require(ggmap)
        require(grid)
        require(maps)
        
        f <- d[d$spp==x,]
        f$prop_max[f$prop_max<10e-4] <- 10e-4 # put a lower limit on the ratios
        f$outlier <- -log10(f$prop_max)
        f <- as.data.frame(f)
        
        # geo space
        coordinates(f) <- c("longitude", "latitude")
        projection(f) <- CRS(prj)
        f <- as.data.frame(spTransform(f, newprj))
        geo <- map + 
                geom_point(data=f[order(f$outlier),], aes(longitude, latitude, color=outlier), size=8) +
                scale_color_gradientn(colours=rev(c("yellow", "orange", "red", "purple", "blue", "black")), 
                                      breaks=seq(0, 3, .5), limits=c(0,3), na.value="white") +
                labs(title=paste0("GEOGRAPHY"),
                     color="maxent\noutlier\nindex") +
                theme(legend.position="left", axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank()) +
                guides(color = guide_colourbar(barheight=20))
        
        # clim space
        fpc <- predict(pc, na.omit(f[,c("cwd", "djf", "jja", "ppt")]))
        pcd <- matrix(ncol=ncol(fpc), nrow=nrow(f))
        pcd[which(!is.na(f$cwd)),] <- fpc
        colnames(pcd) <- colnames(fpc)
        f <- cbind(f, pcd)
        clim <- ggplot() +
                geom_point(data=f[order(f$outlier),], aes(PC1, PC2, color=outlier), size=8) +
                theme_minimal() +
                scale_color_gradientn(colours=rev(c("yellow", "orange", "red", "purple", "blue", "black")), 
                                      breaks=seq(0, 3, .5), limits=c(0,3), na.value="white") +
                coord_fixed() +
                labs(title="CLIMATE SPACE",
                     x="Climate PC1", y="Climate PC2") +
                theme(legend.position="none")
        
        # grid
        grd <- projectRaster(grid, crs=CRS(newprj)) # project grid into lat-long (for plotting over map)
        coordinates(f) <- c("longitude", "latitude")
        projection(f) <- CRS(newprj) 
        grd <- rasterize(f, grd, field="id", fun="count")
        
        txtsz <- theme(text=element_text(size=20))
        
        require(gridGraphics)
        dev.new(width=5, height=5, noRStudioGD=T)
        plot(grd, col=viridis::viridis_pal()(50), axes=F, main="50K DENSITY")
        map("state", add=T)
        grid.echo()
        grid.grab() -> cells
        dev.off()
        p <- arrangeGrob(geo+txtsz, clim+txtsz, rectGrob(gp=gpar(col=NA)), nrow=1)
        p <- arrangeGrob(textGrob(x, gp=gpar(cex=5, fontface="italic")), p, ncol=1, heights=c(.1,1))
        png(paste0(spp_dirs[grepl(x, spp_dirs)], "/outlier_map_v3.png"), width=2400, height=800)
        grid.draw(p)
        pushViewport(viewport(x = .835, y = .44, height = 1, width = .33))    
        grid.draw(cells)
        dev.off()
        
}

# generate plots in parallel
library(doParallel)
#cl <- makeCluster(7)
#registerDoParallel(cl)
results <- foreach(sp = unique(d$spp)) %do% { try(outlier_map(sp)) }
#stopCluster(cl)



### copy charts to second consolidated location

infiles <- list.files("E:/Phylo_modelling/Output/Maxent/V2", pattern="outlier_map_v3.png", recursive=T, full.names=T)

outdir <- "E:/Phylo_modelling/Output/Charts/Mapped_outliers"

for(file in infiles){
      dir <- dirname(file)
      species <- substr(dirname(file), tail(gregexpr("/", dir)[[1]], 1)+1, 1000)
      newfile <- paste0(outdir, "/", species, ".png")
      if(file.exists(newfile)) next()
      #writeLines(paste("moved", species))
      file.copy(file, newfile)
}

