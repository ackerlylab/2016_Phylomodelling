# mapping grid cells for California
## lat long
require(sp)
require(maptools)
require(rgeos)
require(rgdal)
#ca <- readShapeSpatial('/Users/Shared/gisdata/maps/california/CApoly/CA')
ca <- readShapeSpatial('C:/Naia/Data/Background_layers/PROCESSED/GADM_California.shp')
projection(ca) <- CRS('+proj=longlat +datum=WGS84')
plot(ca,axes=T)

cax <- extent(ca)
cax

hdll <- gridlines(ca,easts=seq(-124.5,-114,by=0.5),norths=seq(32,42.5,by=0.5))
plot(ca,col='red')
plot(hdll,add=T)

ca.aea <- spTransform(ca,CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'))
plot(ca.aea)
cax.aea <- extent(ca.aea)
cax.aea

# raster with grid resolution 50x50 km
ras <- raster(xmn=-400000,xmx=550000,ymn=-650000,ymx=450000,res=50000,vals=0,crs=projection(ca.aea))
vals <- cellsFromExtent(ras,ras)
values(ras) <- vals

plot(ras)
plot(ca.aea,add=T)
plot(hd.aea,add=T,col='red')

hdll2aea <- spTransform(hdll,CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'))
plot(ca.aea,col='red')
plot(hdll2aea,add=T)

east_lines <- seq(-350000,500000,by=50000)
north_lines <- seq(-600000,450000,by=50000)
hd.aea <- gridlines(ca.aea,easts=east_lines,norths=north_lines)
class(hd.aea)
plot(ca.aea,col='red',axes=T)
plot(hd.aea,add=T)

par(mfrow=c(1,2))
plot(ca.aea,col='red')
plot(hdll2aea,add=T)

plot(ca.aea,col='red')
plot(hd.aea,add=T)
plot(hdll2aea,add=T)

#dem <- raster('/Users/Shared/gisdata/climate_layers/BCM/CA_topography/ca_270m_t6.asc')
dem <- readRDS('E:/BCM/CA_2014/Summary/HST/Normals_30years/BCM2014_cwd1951-1980_wy_ave_HST.Rdata')[[1]]
dem <- crop(dem,ca.aea)

fx <- c(-4.5e5,6e5)
fy <- c(-7e5,5.5e5)
frame <- SpatialPoints(cbind(fx,fy))

par(mfrow=c(1,1))
#png('/Users/david/Documents/Projects/CalDiversity/PhylodiversityProposal-2012/gridcells/CA_Albers_grids.png',1000,1000)
plot(frame,pch='.',axes=T)
plot(dem,add=T,legend=F)
plot(hd.aea,add=T,col='red')
plot(ca.aea,add=T,lwd=2)
text(5.7e5,north_lines,labels=north_lines/10000)
text(east_lines,4.8e5,labels=east_lines/10000)
text((c(east_lines,max(east_lines)+50000))-25000,-6.3e5,labels=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S'))
text(-4e5,north_lines[1:21]+25000,labels=1:21)
dev.off()

# export lines to kml
setwd('/Users/david/Documents/Projects/CalDiversity/PhylodiversityProposal-2012/gridcells/')
hd.aea2ll <- spTransform(hd.aea,CRS('+proj=longlat +datum=WGS84'))
dd <- data.frame(ID=c('NS','EW'))
row.names(dd) <- c('NS','EW')
hd.aea2ll <- SpatialLinesDataFrame(hd.aea2ll,data=dd)
writeOGR(hd.aea2ll, dsn="Albers_grids.kml", layer="hd.aea2ll", driver="KML")

ca.utm <- spTransform(ca,CRS('+proj=utm +zone=10'))
plot(ca.utm)
(caux <- extent(ca.utm))

## BIODIVERSE GRIDS - 50k
cab <- readShapeSpatial('/Users/david/Documents/Projects/CalDiversity/PhylodiversityProposal-2012/CA_CH_gridcells/test_california/test_california')
plot(cab)
abline(v=0)
abline(h=0)

plot(cab,add=T,col='blue')
plot(ca.aea,add=T)

pl <- SpatialLines2PolySet(hd.aea)
pg <- PolySet2SpatialPolygons(pl)

