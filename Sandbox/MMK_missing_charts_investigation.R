

missing_charts <- c("Artemisia_nesiotica",
                    "Astragalus_ertterae",
                    "Astragalus_nevinii",
                    "Astragalus_traskiae",
                    "Triteleia_clementina",
                    "Brodiaea_kinkiensis",
                    "Calochortus_monanthus",
                    "Camissoniopsis_guadalupensis",
                    "Castilleja_grisea",
                    "Constancea_nevinii",
                    "Cryptantha_traskiae",
                    "Deinandra_clementina",
                    "Dudleya_nesiotica",
                    "Dudleya_traskiae",
                    "Dudleya_virens",
                    "Eriogonum_covilleanum",
                    "Eriogonum_crocatum",
                    "Eriogonum_dasyanthemum",
                    "Eriogonum_davidsonii",
                    "Eschscholzia_ramosa",
                    "Eucephalus_vialis",
                    "Galium_catalinense",
                    "Gambelia_speciosa",
                    "Gilia_nevinii",
                    "Lasthenia_maritima",
                    "Lithophragma_maximum",
                    "Lomatium_insulare",
                    "Lupinus_guadalupensis",
                    "Lupinus_uncialis",
                    "Lycium_verrucosum",
                    "Malacothrix_foliosa",
                    "Malacothrix_junakii",
                    "Munzothamnus_blairii",
                    "Ornithostaphylos_oppositifolia",
                    "Penstemon_barnebyi",
                    "Phacelia_floribunda",
                    "Phacelia_lyonii",
                    "Phlox_muscoides",
                    "Salicornia_rubra",
                    "Senecio_lyonii",
                    "Sibara_filifolia",
                    "Sphaeralcea_munroana",
                    "Trifolium_palmeri",
                    "Zeltnera_davyi",
                    "Zeltnera_trichantha",
                    "Zeltnera_venusta")

missing_charts <- sub("_", " ", missing_charts)
spm <- data.frame(spp=missing_charts, stringsAsFactors=F)

# which are actually missing?
charts <- sub(".png", "", (list.files("E:/Phylo_modelling/Output/Charts/Mapped_outliers")))
spm$missing_charts <- !spm$spp %in% charts
missing_charts <- spm$spp[spm$missing_charts]



##################### run the setup parts of 02_outlier_distributions.r before running this #####################

spm$missing_models <- !spm$spp %in% d$spp
spm$missing_model_inputs <- !spm$spp %in% allSpecies
write.csv(spm, "E:/Phylo_modelling/Output/missing_charts.csv", row.names=F)
spm <- read.csv("E:/Phylo_modelling/Output/missing_charts.csv", stringsAsFactors=F)
maxent_failed <- spm$spp[spm$missing_models==T & spm$missing_model_inputs==F]

stop()




missing_during_charting <- spm$spp[spm$missing_charts==T & spm$missing_models==F]

d_missing <- d[d$spp %in% missing_during_charting,]
table(d_missing$spp)

# point maps of maxent outliers in geographic and climate space
outlier_map <- function(x){
  
  try(dev.off()) 
  
  require(ggplot2)
  require(gridExtra)
  require(raster)
  require(ggmap)
  require(grid)
  
  f <- d[d$spp==x,]
  f$prop_max[f$prop_max<10e-4] <- 10e-4 # put a lower limit on the ratios
  f$outlier <- -log10(f$prop_max)
  f <- as.data.frame(f)
  
  # geo space
  coordinates(f) <- c("longitude", "latitude")
  projection(f) <- CRS(prj)
  f <- as.data.frame(spTransform(f, newprj)) # convert from albers to latlong
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
  grd[1,1] <- 0
  grdd <- as.data.frame(rasterToPoints(grd))
  
  palette <- viridis::scale_fill_viridis(trans="log10", breaks=c(1, 3, 10, 30, 100, 300))
  if(nrow(grdd)==1) palette <- NULL
  
  cells <- map + 
    geom_tile(data=grdd, aes(x, y, fill=layer), color="white") +  
    palette +
    theme_minimal() +
    coord_cartesian() +
    coord_fixed(ratio=1.3) +
    theme(legend.position="right", axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank()) +
    guides(fill = guide_colourbar(barheight=20)) +
    labs(fill="obs/\npx", title="50K DENSITY")
  
  txtsz <- theme(text=element_text(size=20))
  
  # save multiplot
  p <- arrangeGrob(geo+txtsz, clim+txtsz, cells+txtsz, nrow=1)
  p <- arrangeGrob(textGrob(x, gp=gpar(cex=5, fontface="italic")), p, ncol=1, heights=c(.1,1))
  png(paste0(spp_dirs[grepl(x, spp_dirs)], "/outlier_map_v2.png"), width=2400, height=800)
  grid.draw(p)
  dev.off()
  
}

for(x in missing_during_charting){
  print(x)
  try(outlier_map(x))
}


######################## run 01_simple_maxent_modeling setup before running this chunk ##############################




for(i in 1:length(missing_charts)) {
  options(java.parameters = "-Xmx1g" )
  require(dismo)
  require(rJava)
  require(raster) 
  
  mySpecies <- missing_charts[i]
  print(mySpecies)
  #Sys.time()
  
  # Read in species occurrence data and background
  pres <- readRDS(paste(spdir, mySpecies, ".rdata", sep=""))
  
  if(nrow(pres)==0){ ##################################################
    print("           0 records in input file")
    next()
  }
  
  coordinates(pres) <- ~longitude + latitude
  projection(pres) <- CRS(orig.project)
  occur <- spTransform(pres, projection(predictors))
  
  #extract( predictors[[1]], occur)
  
  # Directory to write files to
  #mx.dir = paste(odir, 'Maxent/V1/',mySpecies, sep="")
  mx.dir = paste(odir, 'Maxent/V2/',mySpecies, sep="")
  if(file.exists(mx.dir) == F) {dir.create(mx.dir, recursive=T)} 
  
  occur_raw <- occur
  occur <- occur[!is.na(extract(predictors[[1]], occur)),] ###############################################
  if(nrow(occur)==0){
    print("           no occurrences overlapping climate data")
    next()
  }
  if(nrow(occur)/nrow(occur_raw) < .5) print("           FIXED: more than half of occurences had climate NAs")
  
  # Run the model!
  mx <- try(maxent(predictors, p=occur, a=bg ,progress='text', path=mx.dir, args=mxArgs))
  saveRDS(mx, file = paste(mx.dir, 'ModelObject.rdata', sep='/'))
  
  
  # Predict suitability at occurrences
  cp <- extract(predictors, occur)
  
  predp <- predict(mx, cp)
  preda <- predict(mx, mx@absence)
  
  # Save results to .rdata object
  results <- list()
  results[['summary']] <- data.frame(species=mySpecies, auc.train=as.numeric(mx@results['Training.AUC',]),records.used=nrow(mx@presence))
  results[['occur']] <- occur
  results[['occur']]$suitability <- predp
  results[['occur']] <- cbind(results$occur, as.data.frame(cp))
  results[['summary']]$max.suit.occ <- max(predp, na.rm=T)
  results[['bg']] <- bg
  results[['bg']]$bgsuitability <- preda
  results[['summary']]$max.suit.all <- max(c(predp, preda), na.rm=T)
  results[['occur']]$prop_max.suit.occ <- results[['occur']]$suitability / results[['summary']]$max.suit.occ
  results[['occur']]$prop_max.suit.all <- results[['occur']]$suitability / results[['summary']]$max.suit.all
  
  saveRDS(results, paste(mx.dir, 'Suitability_data.rdata',sep='/'))
}







