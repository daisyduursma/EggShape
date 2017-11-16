
library(raster)
#library(RColorBrewer)
#library(lattice)
library(maptools)
library(letsR)
#library(plantecophys)
library(data.table)

traits<-read.csv('/Users/daisy/Google Drive/PhD/ThermalStress/data/ThermalStressTraits2016-05-09.csv')
traits$BL.name<-with(traits,paste0("/",BLGenus,"_",BLSpecies,"_",number))
traits<-subset(traits,USEThermalStressStudy != 'omit')

passer<-traits$Species

obs <- fread("~/Google Drive/PhD/Data/Observaitons/Cleaned/observations/PASSARINE_ALAGBIFABBBSATLASGridCellIDObservations2016-06-17.csv",data.table=FALSE)

gnt<-subset(read.csv('~/Google Drive/PhD/birdTraits/PublishedData/Australian_Bird_Data_Version_1.csv'),
            is.na(X6_Subspecies_name_2))

gnt <- subset(read.csv('~/Google Drive/PhD/birdTraits/PublishedData/Australian_Bird_Data_Version_1.csv'),
              X31_Non.breeding_populations_of_core_taxa_4==0 &
                X32_Extinct_4==0 &
                X34_Vagrant_4==0 &
                X17_Species_2==1 &
                is.na(X6_Subspecies_name_2) &
                X11_Order_2=="Passeriformes")[,c("Species","X193_National_movement_local_dispersal_13")]


############################

#100 km shapefile mask based on climate data using in study

############################
aus<-raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/BOM/maxann.txt")#mask of Australia
aus_XY <- as.data.frame(rasterToPoints(aus))#find lat and long of 1km data
aus<-aggregate(aus,fact=40,na.rm=TRUE)#make raster 1 degree
aus[]<-1:1505 #Give each grid cell unique value
aus_XY$GridCellID<-extract(aus,cbind(aus_XY$x,aus_XY$y),method='simple')#extract the 1km resolution cell value
gridCount<-as.data.frame(table(aus_XY$GridCellID))#find how many grid cells fall in the 1 degree data
gridCount<-subset(gridCount,Freq>=533)#only keep grid cells that are actually at least 33% land
gridCount$val<-1 #set the value to 1
ausMask<-subs(aus, data.frame(gridCount[,c("Var1","val")])) 
aus<-raster::mask(aus,ausMask)



#locations of shapefiles
shape.dir<-"~/Google Drive/PhD/Data/Spatial/BirdlifeSpeciesDistributions/All"


#function to make polygon 
Breedingpoly<-function(shape){
  #Breeding polygon
  TF <- (shape$SEASONAL == 2)
  if (length(shape$SEASONAL) >= 2) {#i<-51 has both breeding and resident 
    if (length(TF[TF == TRUE]) >= 1) {
      breeding <-
        lets.shFilter(shape, presence = "1", origin = "1", seasonal = "2")
      breeding <- crop(shape,extent(aus))
    }
  }
  #year round resident population
  TFres <- (shape$SEASONAL == 1)
  if (length(TFres[TFres == TRUE]) >= 1) {
    resident <-
      lets.shFilter(shape, presence = "1", origin = "1", seasonal = "1")
    resident <- crop(shape,extent(aus))
  }
  #combine polygons to one shapefile, works even if only breeding or resident polygon
  if(length(TFres[TFres == TRUE]) >= 1 & length(TF[TF == TRUE]) >= 1){
    Bpoly<-merge(breeding,resident,by = intersect(names(x), names(y)))
  } else{
    if(length(TFres[TFres == TRUE]) >= 1){#just resident
      Bpoly<-resident 
    } else {
      if (length(TF[TF == TRUE]) >= 1){#just breeding
        Bpoly<-breeding
      }
    }
  }
  return(Bpoly)
}

#Function to find unique 1 degree grid cells with presence
griddedObs<-function(spobsXY,ausMask) {
  obsXYPoints <-
    SpatialPoints(
      cbind(spobsXY$lon,spobsXY$lat),
      proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'))
  obsGridcell <-
    na.omit(as.data.frame(extract(ausMask,obsXYPoints)))#find the 100km cells occurs in (based on centerpoint of grid cells) for observations
  obsGridcell$value = 1
  colnames(obsGridcell) <- c("ID","value")
  obsGridcell <- unique(obsGridcell)
  obsRaster <- subs(ausMask,
                    data.frame(obsGridcell[,c("ID","value")]))
  specObsXY <- as.data.frame(rasterToPoints(obsRaster))
  return (specObsXY)
}

#species that I use the points from 
spPoints<-subset(traits,USEThermalStressStudy=="points")
spClip<-subset(traits,USEThermalStressStudy=="clip")
spManual<-subset(traits,USEThermalStressStudy!="clip"&USEThermalStressStudy!="points")

###########################

#Species using points with

##########################
obsGridPoints<-list()#list of observations
obsPoints<-list()#list of gridded observations
for(i in 1:nrow(spPoints)) {
  #prep observaitons
  spobs <-subset(obs,Scientific.Name == as.vector(spPoints$Species[i]))# get species observation
  movement<-subset(gnt,Species==as.character(spPoints$Species[i]), 
                   X193_National_movement_local_dispersal_13)#find out if local dispersion, if movement is not limited to local then use only breeding observations
  if(movement$X193_National_movement_local_dispersal_13!=1){ 
    spobs<-subset(spobs,Breeding==1)
  }
  #convert the data to presence in gridcell
  spGrid<-griddedObs(spobs,aus)
  spGrid$Species<-as.vector(spPoints$Species[i])
  obsPoints[[i]]<-spobs
  obsGridPoints[[i]]<-spGrid
}
#turn the data back into a dataframe
obsPoints<-do.call("rbind",obsPoints)
obsGridPoints<-do.call("rbind",obsGridPoints)

  
###########################

#Species clipping to polygon

##########################
  
obsGridClip<-list()
obsClip<-list()
for(ii in 1:nrow(spClip)) {
  #prep observaitons
  spobsClip <-subset(obs,Scientific.Name == as.vector(spClip$Species[ii]))# get species observation
  movement<-subset(gnt,Species==as.character(spClip$Species[ii]), 
                   X193_National_movement_local_dispersal_13)#find out if local dispersion, if movement is not limited to local then use only breeding observations
  if(movement$X193_National_movement_local_dispersal_13!=1){ 
    spobsClip<-subset(spobsClip,Breeding==1)
  }
  #clip observation to birdlife Polygon
  Polyname<-subset(traits,Species==as.vector(spClip$Species[ii]),BL.name)
  s <- readShapePoly(paste0(shape.dir,#read in species shape
                            Polyname$BL.name),
                     proj4string = CRS("+proj=longlat +datum=WGS84"))
  Bpoly<-Breedingpoly(s)
  specR <- rasterize(Bpoly,aus,getCover = TRUE)#rasterize shapefile
  cfun <-
    function(x) {
      x[x < 1] <- NA; return(x)
    }#set background to NA
  specR <- calc(specR, cfun)
  specR <- (specR / specR) * aus#clip to australia boundaries
  s2 <-
    buffer(specR,300000,dissolve = TRUE) * aus / aus#buffer around shapefile to further clean polygons
  #identify obs inside buffered polygon
  spobsClip$insidePoly <- #identify obs inside buffered polygon
    as.vector(extract(s2,cbind(spobsClip$lon,spobsClip$lat)))
  spobsClip <- subset(spobsClip,insidePoly == 1)#inside
  #convert the data to presence in gridcell
  spGridfromClip<-griddedObs(spobsClip,aus)
  spGridfromClip$Species<-as.vector(spClip$Species[ii])
  obsClip[[ii]]<-spobsClip
  obsGridClip[[ii]]<-spGridfromClip
  message(ii)
  if(nrow(spGridfromClip == 0) == TRUE){ stop("no gridded observations")}
}
#turn the data back into a dataframe
obsClip<-do.call("rbind",obsClip)
obsGridClip<-do.call("rbind",obsGridClip)


###########################

#Species with manual clipping of observations

##########################


spManual<-spManual[,c("Species","USEThermalStressStudy")]

# Heteromyias albispecularis                   points between -25and -10
subset(gnt,Species=="Heteromyias albispecularis", 
                 X193_National_movement_local_dispersal_13)

HA<-subset(obs,Scientific.Name == "Heteromyias albispecularis" & 
             lat >=-25 &
             lat <=-10)
plot(aus)
plot(SpatialPoints(
  cbind(HA$lon,HA$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)



# 52            Myiagra inquieta               points but need to remove one
subset(gnt,Species=="Myiagra inquieta", 
       X193_National_movement_local_dispersal_13)

MI<-subset(obs,Scientific.Name == "Myiagra inquieta" 
           & lat >= -40)
plot(aus)
plot(SpatialPoints(
    cbind(MI$lon,MI$lat),
    proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)




# 54          Myiagra ruficollis              everything between -25 and -10

subset(gnt,Species=="Myiagra ruficollis", 
       X193_National_movement_local_dispersal_13)

MR<-subset(obs,Scientific.Name == "Myiagra ruficollis" & 
             lat >=-25 &
             lat <=-10)
plot(aus)
plot(SpatialPoints(
  cbind(MR$lon,MR$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)

# 66         Petroica multicolor          points, western population correct
subset(gnt,Species=="Petroica multicolor", 
       X193_National_movement_local_dispersal_13)

PM<-subset(obs,Scientific.Name == "Petroica multicolor" & 
             lon <=125)
plot(aus)
plot(SpatialPoints(
  cbind(PM$lon,PM$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)



# 67         Philemon buceroides                     everything north of -25
subset(gnt,Species=="Philemon buceroides", 
       X193_National_movement_local_dispersal_13)

PB<-subset(obs,Scientific.Name == "Philemon buceroides" & 
             lat >=-25)
plot(aus)
plot(SpatialPoints(
  cbind(PB$lon,PB$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)



# 82         Tregellasia leucops             points north of -20 through -10
subset(gnt,Species=="Tregellasia leucops", 
       X193_National_movement_local_dispersal_13)

TL<-subset(obs,Scientific.Name == "Tregellasia leucops" & 
             lat >=-20 &
             lat <=-10)
plot(aus)
plot(SpatialPoints(
  cbind(TL$lon,TL$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)

# 86       Zosterops citrinellus             points north of -20 through -10

subset(gnt,Species=="Zosterops citrinellus", 
       X193_National_movement_local_dispersal_13)

ZC<-subset(obs,Scientific.Name == "Zosterops citrinellus" & 
             lat >=-20 &
             lat <=-10)
plot(aus)
plot(SpatialPoints(
  cbind(ZC$lon,ZC$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)




# 96          Acanthiza apicalis                      points -40 through -10
subset(gnt,Species=="Acanthiza apicalis", 
       X193_National_movement_local_dispersal_13)

AA<-subset(obs,Scientific.Name == "Acanthiza apicalis" & 
             lat >=-40 &
             lat <=-10)
plot(aus)
plot(SpatialPoints(
  cbind(AA$lon,AA$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)


# 98           Acanthiza ewingii                     points only in tasmania
subset(gnt,Species=="Acanthiza ewingii", 
       X193_National_movement_local_dispersal_13)

AE<-subset(obs,Scientific.Name == "Acanthiza ewingii" & 
             lat <=-40) 
plot(aus)
plot(SpatialPoints(
  cbind(AE$lon,AE$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)


# 124    Anthochaera chrysoptera                  western population correct
subset(gnt,Species=="Anthochaera chrysoptera", 
       X193_National_movement_local_dispersal_13)

ACHa<-subset(obs,Scientific.Name == "Anthochaera chrysoptera" & 
             lon <=119)
ACHb<-subset(obs,Scientific.Name == "Anthochaera chrysoptera" & 
               lon >=135)
ACH<-rbind(ACHa,ACHb)
ACH<-subset(ACH, 
        lat <=-20)

plot(aus)
plot(SpatialPoints(
  cbind(ACH$lon,ACH$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)


# 138       Atrichornis clamosus                            points 110 - 120

subset(gnt,Species=="Atrichornis clamosus", 
       X193_National_movement_local_dispersal_13)

AC<-subset(obs,Scientific.Name == "Atrichornis clamosus" & 
             lon >=110 &
             lon <=120)
plot(aus)
plot(SpatialPoints(
  cbind(AC$lon,AC$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)


# 241     Menura novaehollandiae points, 140 through 160 and -25 through -50
subset(gnt,Species=="Menura novaehollandiae", 
       X193_National_movement_local_dispersal_13)

MN<-subset(obs,Scientific.Name == "Menura novaehollandiae" & 
             lon >=140 &
             lat <=-25)
plot(aus)
plot(SpatialPoints(
  cbind(MN$lon,MN$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)


# 244         Neochmia ruficauda                      points -40 through -10
subset(gnt,Species=="Neochmia ruficauda", 
       X193_National_movement_local_dispersal_13)

NR<-subset(obs,Scientific.Name == "Neochmia ruficauda" & 
             lat >=-40 &
             lat <=-10)
plot(aus)
plot(SpatialPoints(
  cbind(NR$lon,NR$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)


# 260         Petroica phoenicea             points -25 to -50 and 135 - 160
subset(gnt,Species=="Petroica phoenicea", 
       X193_National_movement_local_dispersal_13)

PP<-subset(obs,Scientific.Name == "Petroica phoenicea" & 
             lon >=135 &
             lat <=-25)
plot(aus)
plot(SpatialPoints(
  cbind(PP$lon,PP$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)

# 272  Poecilodryas superciliosa                       points -30 though -10

subset(gnt,Species=="Poecilodryas superciliosa", 
       X193_National_movement_local_dispersal_13)

PS<-subset(obs,Scientific.Name == "Poecilodryas superciliosa" & 
             lat >=-30 &
             lat <=-10)
plot(aus)
plot(SpatialPoints(
  cbind(PS$lon,PS$lat),
  proj4string = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')),add=TRUE)




obsManual<-rbind(HA,MI,MR,PM,PB,AA,ZC,TL,AE,ACH,AC,MN,NR,PP,PS)

obsGridManual<-list()
for(j in 1:nrow(spManual)) {
  #prep observaitons
  spobs <-subset(obsManual,Scientific.Name == as.vector(spManual$Species[j]))# get species observation
  movement<-subset(gnt,Species==as.vector(spManual$Species[j]), 
                   X193_National_movement_local_dispersal_13)#find out if local dispersion, if movement is not limited to local then use only breeding observations
  if(movement$X193_National_movement_local_dispersal_13!=1){ 
    spobs<-subset(spobs,Breeding==1)
  }
  #convert the data to presence in gridcell
  spGrid<-griddedObs(spobs,aus)
  spGrid$Species<-as.vector(spManual$Species[j])
  obsGridManual[[j]]<-spGrid
}
#turn the data back into a dataframe
obsGridManual<-do.call("rbind",obsGridManual)

cleanObs<-unique(rbind(obsClip[,c("Scientific.Name","lat","lon", "Breeding")],obsPoints,obsManual))

cleangriddedObs<-unique(rbind(obsGridManual,obsGridClip,obsGridPoints))

write.csv(cleanObs, paste0('~/Google Drive/PhD/ThermalStress/data/PassarineObservationsFinal',
                                as.Date(Sys.time()),'.csv'),row.names=FALSE)

write.csv(cleangriddedObs, paste0('~/Google Drive/PhD/ThermalStress/data/PassarineObservationsGriddedFinal',
                           as.Date(Sys.time()),'.csv'),row.names=FALSE)






