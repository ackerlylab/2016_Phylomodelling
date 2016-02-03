# Script for running MaxEnt models for California plants

# Clear workspace
rm(list=ls())

#-----------------#
# Set directories #
#-----------------#

# source user_parameters.r before running

wdir <- project_stem_dir
cdir <- climate_data_dir

spdir <- paste(wdir, '/Data/Species/', sep='')
outspdir <- paste(wdir, 'Data/Species/Processed3/',sep='/')

#----------------#
# Load libraries #
#----------------#
require(raster)

#----------------#
# Set parameters #
#----------------#
#what columns of interest in input species csv files?
#cols <- c('id','current_name_binomial', 'latitude', 'longitude') #for data filtering step
cols <- c('id','current_name_binomial', 'latitude', 'longitude', 'datum', 'exclude_species') #for using cleaned data in final models
#what climate files and variables for modelling?
climnames <- sort(c('cwd','djf','jja','ppt'))
env.files <- list.files(path=cdir, pattern='data', full.names=FALSE)


#------------------------#
# Prepare climate layers #
#------------------------#
# read in raw climate layers
files <- env.files[which(substr(env.files,9,11)%in%climnames)]
predictors <- stack(lapply(files,function(x) readRDS(paste(cdir,x,sep=""))))
names(predictors) = climnames 

# # fill in gaps in climate layers (~1 km coast line for all and also streams and lakes for cwd - only for flagging!)
# for(a in 2:3) {
#   lyr <- predictors[[a]]
#   lyr_fill <- focal(lyr, w=matrix(1,nrow=3,ncol=3), function(x)mean(x[!is.na(x)]), NAonly=T)
#   for(i in 1:4) lyr_fill <- focal(lyr_fill, w=matrix(1,nrow=3,ncol=3), function(x)mean(x[!is.na(x)]), NAonly=T)
#   saveRDS(lyr_fill,paste(cdir,climnames[a],'_filled.rdata',sep=''))
#   print(Sys.time())
# }

# read in filled climate layers
files2 <- env.files[which(substr(env.files,5,10)%in%'filled')]
fill_predictors <- stack(lapply(files2,function(x) readRDS(paste(cdir,x,sep=""))))
names(fill_predictors) = climnames


# Aggregate climate layers to ~1km resolution
# for(i in 1:4) {
#   lyr <- predictors[[i]]
#   ag <- aggregate(lyr,4,na.rm=T,fun=mean)
#   saveRDS(ag, paste(cdir, climnames[i],'_1080km.rdata', sep=''))
# }
#preds1080 <- stack(lapply(climnames,function(x) readRDS(paste(cdir, x,'_1080m.rdata', sep=''))))

# Aggregate climate layers to ~800m resolution
for(i in 1:4) {
      lyr <- predictors[[i]]
      ag <- aggregate(lyr,3,na.rm=T,fun=mean)
      saveRDS(ag, paste(cdir, climnames[i],'_810m.rdata', sep=''))
}
preds810 <- stack(lapply(climnames,function(x) readRDS(paste(cdir, x,'_810m.rdata', sep=''))))


#---------------------------#
# Prepare 50x50 km aea grid #
#---------------------------#
ca <- readShapeSpatial('C:/Naia/Data/Background_layers/PROCESSED/GADM_California.shp')
projection(ca) <- CRS('+proj=longlat +datum=WGS84')
ca.aea <- spTransform(ca,CRS('+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'))
cax.aea <- extent(ca.aea)
cax.aea

# raster with grid resolution 50x50 km
ras <- raster(xmn=-400000,xmx=550000,ymn=-650000,ymx=450000,res=50000,vals=0,crs=projection(ca.aea))

#east_lines <- seq(-350000,500000,by=50000)
#north_lines <- seq(-600000,450000,by=50000)
#hd.aea <- gridlines(ca.aea,easts=east_lines,norths=north_lines)
#class(hd.aea)
#plot(ras,axes=T)
#plot(ca.aea,add=T)
#plot(hd.aea,col='red',add=T)

#saveRDS(ras, paste(cdir, 'background50by50km.Rdata', sep=''))

#---------------------------#
# Prepare occurrence tables #
#---------------------------#
files <- list.files(path=spdir, pattern='.csv', full.names=T)
df <- do.call("rbind", lapply(files, FUN=function(x){
      f <- read.csv(x,as.is=T)[,cols]
      f$source_file <- x
      return(f)
}))

splist_raw <- sort(unique(df[,'current_name_binomial']))

# remove 'species' with void name or non-binomial name
splist <- splist_raw[grep(' ',splist_raw)]
#saveRDS(splist, file=paste(outspdir,'0_Species_list.rdata',sep=''))

# clean occurrences
df_clean <- subset(df, current_name_binomial%in%splist & !is.na(latitude) & !is.na(longitude) & latitude >0)
row.names(df_clean) <- NULL
splist_clean <- sort(unique(df_clean$current_name_binomial))
saveRDS(splist_clean, file=paste(outspdir,'0_Species_list_v2.rdata',sep=''))

#For each species, save dataframe with occurrences
for(i in 1:length(splist_clean)) {
      species <- splist_clean[i]
      occur_raw <- subset(df_clean, current_name_binomial==species)
      occur <- subset(occur_raw, !is.na(latitude) & !is.na(longitude))
      row.names(occur) <- NULL
      
      saveRDS(occur, file=paste(outspdir,species,'.rdata',sep=''))
      print(i)
}


#----------------------------#
# How many records in water? #
#----------------------------#
# Get tabulated data % records for each species that are in the water and
# how many are preserved by using filled climate layers

# project lat long to same projection as climate layers
orig.project <- '+proj=longlat +ellps=WGS84'
pts <- df_clean
coordinates(pts) <- ~longitude + latitude
projection(pts) <- CRS(orig.project)
p.pts <- spTransform(pts, projection(predictors))

df_clean$hasClim <- 1
vals <- extract(predictors[[1]], p.pts)
df_clean$hasClim[is.na(vals)] <- 0

df_clean$hasClimFill <- 1
vals2 <- extract(fill_predictors[[1]], p.pts)
df_clean$hasClimFill[is.na(vals2)] <- 0

df_clean$hasClim1080 <-1
vals3 <- extract(preds1080[[1]], p.pts)
df_clean$hasClim1080[is.na(vals3)] <- 0

tab <- data.frame(splist=sort(unique(df_clean$current_name_binomial)))
tab$splist <- as.character(tab$splist)
tab$records <- as.numeric(table(df_clean$current_name_binomial))
tab$hasClim[match(names(table(df_clean[df_clean$hasClim==1,'current_name_binomial'])),tab$splist)] <- as.numeric(table(df_clean[df_clean$hasClim==1,'current_name_binomial']))


tab$hasClimFill[match(names(table(df_clean[df_clean$hasClimFill==1,'current_name_binomial'])),tab$splist)] <- as.numeric(table(df_clean[df_clean$hasClimFill==1,'current_name_binomial']))

tab$hasClim1080[match(names(table(df_clean[df_clean$hasClim1080==1,'current_name_binomial'])),tab$splist)] <- as.numeric(table(df_clean[df_clean$hasClim1080==1,'current_name_binomial']))
#replace NA with 0
tab$hasClim[is.na(tab$hasClim)] <- 0
tab$hasClimFill[is.na(tab$hasClimFill)] <- 0
tab$hasClim1080[is.na(tab$hasClim1080)] <- 0

tab$records <- as.numeric(table(df_clean$current_name_binomial))

#saveRDS(tab,paste(outspdir, '00_Tabulation_records_with_climate.rdata',sep=''))



# Plotting
# plot(djf)
# points(p.pts[is.na(vals2),],col='blue')
# df_final2 <- df_clean[!is.na(vals),]
# 
# #ColRiv=drawExtent()
# PR=drawExtent()
# plot(crop(cwd,PR))
# points(p.pts[is.na(vals),],col='blue')


#--------------------------#
# Prepare background table #
#--------------------------#
# Get climate variables for each occurrence record

# clean occurrences
#df_clean <- subset(df, current_name_binomial%in%splist & !is.na(latitude) & !is.na(longitude) & latitude >0)
#row.names(df_clean) <- NULL

# read in climate layers
files <- env.files[which(substr(env.files,9,11)%in%climnames)]
#predictors <- stack(lapply(files,function(x) readRDS(paste(cdir,x,sep=""))))
#predictors <- stack(lapply(climnames, function(x) readRDS(paste(cdir, x, '_1080m.rdata',sep=''))))
predictors <- stack(lapply(climnames, function(x) readRDS(paste(cdir, x, '_810m.rdata',sep=''))))
names(predictors) = climnames 

# project lat long to same projection as climate layers
orig.project <- '+proj=longlat +ellps=WGS84'
pts <- df_clean
coordinates(pts) <- ~longitude + latitude
projection(pts) <- CRS(orig.project)
p.pts <- spTransform(pts, projection(predictors))

# get unique cell ids with records
cells <- unique(cellFromXY(predictors[[1]],p.pts))
points <- xyFromCell(predictors[[1]],cells,spatial=T)
vals <- extract(predictors[[1]],points)

notna <- points[!is.na(vals)]

# sample 10,000 cells to be used for background
bg <- sample(notna,10000)
#saveRDS(bg, file=paste(outspdir,'10000_CA_plants_bg.rdata',sep=''))
#saveRDS(bg, file=paste(outspdir,'10000_CA_plants_bg_1080m.rdata',sep=''))
saveRDS(bg, file=paste(outspdir,'10000_CA_plants_bg_810m.rdata',sep=''))


########## alternative method: sample background points to match occurrence locality distribution #########

# load climate data (used to eliminate water occurrences)
climfiles <- list.files(path=cdir, pattern='data', full.names=T)
climfiles <- env.files[grepl(paste(paste0(climnames, "_filled"), collapse="|"), env.files)]
library(raster)
climrast <- do.call("stack", lapply(climfiles, readRDS))
climrast <- sum(climrast)

# load occurrences
allocc <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Species/processed3/atomic", full.names=T)
allocc <- lapply(allocc, readRDS)
allocc <- do.call("rbind", allocc)
allocc <- allocc[!is.na(allocc$longitude + allocc$latitude),]
coordinates(allocc) <- c("longitude", "latitude")
projection(allocc) <- '+proj=longlat +ellps=WGS84'
allocc <- spTransform(allocc, crs(climrast))

# sample 10k records for background, export result
bg <- allocc[!is.na(extract(climrast, allocc)),]
set.seed(123)
bg <- landocc[sample(length(bg), 10000),]
saveRDS(bg, file=paste(outspdir,'10000_CA_plants_bg_810m_occbias.rdata',sep=''))



#---------------------------#
# Flag species not modelled #
#---------------------------#
df_bad <- subset(df, !current_name_binomial%in%splist | is.na(latitude) | is.na(longitude) | latitude <0)[,c(
      'id','current_name_binomial','latitude','longitude')]
df_bad$flag <- NA
df_bad$flagreason <- NA
df_bad$flag[is.na(df_bad$latitude) | is.na(df_bad$longitude)] <- 1
df_bad$flagreason[df_bad$flag == 1] <- 'lat or long missing'

df_bad$flag[df_bad$latitude < 0 | df_bad$longitude > 0] <- 2
df_bad$flagreason[df_bad$flag == 2] <- 'negative lat or positive long'

df_bad$flag[!df_bad$current_name_binomial %in% splist] <- 3
df_bad$flagreason[df_bad$flag == 3] <- 'invalid species binomial'

df_bad$flag[df_bad$current_name_binomial == ''] <- 4
df_bad$flagreason[df_bad$flag == 4] <- 'missing species binomial'

df_bad$flag[!df_bad$current_name_binomial %in% splist & (is.na(df_bad$latitude) | df_bad$latitude < 0) |
                  !df_bad$current_name_binomial %in% splist & (is.na(df_bad$longitude) | df_bad$longitude > 0)] <- 5
df_bad$flagreason[df_bad$flag == 5] <- 'multiple errors'

row.names(df_bad) <- NULL
# extra bad occurrences based on climate coverage
df_bad2 <- subset(df_clean,hasClim==0 | hasClim1080==0)
df_bad2$flag <- df_bad2$flagreason <- NA
df_bad2$flag[df_bad2$hasClim==0] <- 6
df_bad2$flagreason[df_bad2$flag == 6] <- 'missing raw climate but modelled'

df_bad2$flag[df_bad2$hasClim1080==0] <- 7
df_bad2$flagreason[df_bad2$flag == 7] <- 'in water or out of state'

df_bad_all <- rbind(df_bad,df_bad2[,c('id','current_name_binomial','latitude','longitude','flag','flagreason')])

write.csv(df_bad_all,paste(outspdir, '00_flagged_occurrences_not_modelled.csv',sep=''))
