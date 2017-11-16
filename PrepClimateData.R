# make layers for Annual, spring, summer VPD
#extract mean, median, and variance for rasterized breeding range of species
#climate to extract: VPD (3 time periods)
#                     annual temperture 
#                     annual wind speeds
#                     EVI
#                     FPAR 
#                     distance from standing water


#############################

#MAKE VPD

#############################
rm(list = ls())
library(raster)

work.dir<-'/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/BOM'
shape.dir<-"~/Google Drive/PhD/Data/Spatial/BirdlifeSpeciesDistributions/All"

traits<-read.csv('/Users/daisy/Google Drive/PhD/ThermalStress/data/ThermalStressTraits2016-05-09.csv')
traits$BL.name<-with(traits,paste0("/",BLGenus,"_",BLSpecies,"_",number))#make shapefile name


#make climate data

rhSpring<-raster(paste0(work.dir,"/rh15sep.txt" ))+
    raster(paste0(work.dir,"/rh15oct.txt" ))+
    raster(paste0(work.dir,"/rh15nov.txt" ))
rhSpring<-aggregate(rhSpring,fact=5,fun=mean)

rhSummer<-raster(paste0(work.dir,"/rh15dec.txt" ))+
  raster(paste0(work.dir,"/rh15jan.txt" ))+
  raster(paste0(work.dir,"/rh15feb.txt" ))
rhSummer<-aggregate(rhSummer,fact=5,fun=mean)

rhann<-raster(paste0(work.dir,"/rh15an.txt" ))
rhann<-aggregate(rhann,fact=5,fun=mean)

tmaxAn<-raster(paste0(work.dir,"/maxann.txt" ))
tmaxAn<-aggregate(tmaxAn,fact=5,fun=mean)

tmaxSpring<-raster(paste0(work.dir,"/maxspr.txt" ))
tmaxSpring<-aggregate(tmaxSpring,fact=5,fun=mean)

tmaxSummer<-raster(paste0(work.dir,"/maxsum.txt" ))
tmaxSummer<-aggregate(tmaxSummer,fact=5,fun=mean)

#make a mask
aus<-raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/BOM/maxann.txt")#mask of Australia
aus_XY <- as.data.frame(rasterToPoints(aus))#find lat and long of 1km data
aus<-aggregate(aus,fact=20,na.rm=TRUE)#make raster .5 degree
aus[]<-1:5865 #Give each grid cell unique value
aus_XY$GridCellID<-extract(aus,cbind(aus_XY$x,aus_XY$y),method='simple')#extract the 1km resolution cell value
gridCount<-as.data.frame(table(aus_XY$GridCellID))#find how many grid cells fall in the 1 degree data
gridCount<-subset(gridCount,Freq>=134)#only keep grid cells that are actually at least 33% land
gridCount$val<-1 #set the value to 1
ausMask<-subs(aus, data.frame(gridCount[,c("Var1","val")])) 
aus<-raster::mask(aus,ausMask)

