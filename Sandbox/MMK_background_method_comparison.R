


# compare effects of maxent background schemes on richness maps

library(raster)
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

for(min in c(0,10,30)){
        a <- raster(paste0("C:/Lab_projects/2016_Phylomodelling/Output/Richness/V4/richness_810m_min", min, ".tif"))
        b <- raster(paste0("C:/Lab_projects/2016_Phylomodelling/Output/Richness/V5/richness_810m_min", min, ".tif"))
        d <- b-a
        
        s <- stack(a, b, d)
        s <- as.data.frame(rasterToPoints(s))
        names(s) <- c("x", "y", "cells", "records", "difference")
        
        for(stat in c("species", "percent")){
                
                title <- "Species richness difference:\n (recordMethod - cellMethod)"
                if(stat=="percent"){
                        s$difference <- s$difference / s$cells
                        title <- "Species richness difference:\n (recordMethod - cellMethod) / cellMethod"
                }
                
                p <- ggplot(s, aes(x, y, fill=difference)) +
                        geom_raster() +
                        scale_fill_gradientn(colours=c("red", "darkred", "black", "darkgreen", "green"),
                                             limits=max(abs(range(na.omit(s$difference)))) * c(-1,1)) +
                        coord_fixed(ratio=1) +
                        theme(panel.background=element_blank(), panel.grid=element_blank(),
                              axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                              legend.position=c(.7, .92), legend.direction="horizontal", title=element_text(size=25)) +
                        guides(fill = guide_colorbar(barwidth=30)) +
                        labs(title=title, fill=NULL)
                
                png(paste0("C:/Lab_projects/2016_Phylomodelling/Output/Richness/V5/richness_difference_n", min, "_v4_v5_", stat, ".png"), width = 1200, height = 1500)
                plot(p)
                dev.off()
        }
}

