#calculate completeness index for 50km gridcells and 100km gridcells
#and remove observations most likely to be erroneous

rm(list = ls())
library(raster)
library(Hmisc)
library(data.table)
#library(rgeos)
library(SpatialTools)
library(RColorBrewer)
library(gtools)
library(vegan)

Albers<-"+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

km1<-raster(nrows=4900, ncols=5600, 
            xmn=-2904860, xmx=2695140,ymn=-5334205, ymx=-434205,
            crs =paste(Albers))



#read in observations
occ <- fread("~/Google Drive/PhD/ThermalStress/data/PassarineObservationsFinal2016-06-21.csv",
                        data.table=FALSE)
species <- unique(occ$Scientific.Name)

#change XY to albers equal area
Locs_sp <- SpatialPoints(cbind(occ$lon, occ$lat), proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
Locs_sp <- spTransform(Locs_sp, CRS=CRS(paste(Albers)))
Alberslocs<-as.data.frame(Locs_sp)
colnames(Alberslocs) <- c("Lon_Albers","Lat_Albers")
occ<-cbind(occ,Alberslocs)



aus<-raster::raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/BOM/maxann.txt",
                    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")#mask of Australia
aus<-aus/aus #make it into a maks
ausAlbers<-projectRaster(aus,km1,method='ngb')#project to albers
aus_XY <- as.data.frame(rasterToPoints(ausAlbers))#find lat and long of 1km data
aus50<-aggregate(ausAlbers,fact=50,na.rm=TRUE)#
aus100<-aggregate(ausAlbers,100,'mean')#100kmX100km resolution
aus50[]<-1:10976
aus100[]<-1:2744 #Give each grid cell unique value
aus_XY$GridCell50ID<-raster::extract(aus50,cbind(aus_XY$x,aus_XY$y),method='simple')#extract the 1km resolution cell value
aus_XY$GridCell100ID<-raster::extract(aus100,cbind(aus_XY$x,aus_XY$y),method='simple')#extract the 1km resolution cell value
gridCount50<-as.data.frame(table(aus_XY$GridCell50ID))#find how many grid cells fall in the 1 degree data
gridCount100<-as.data.frame(table(aus_XY$GridCell100ID))
gridCount50<-subset(gridCount50,Freq>=833)#only keep grid cells that are actually at least 33% land
gridCount50$val<-1 #set the value to 1
gridCount100<-subset(gridCount100,Freq>=3333)#only keep grid cells that are actually at least 33% land
gridCount100$val<-1 #set the value to 1
aus100Mask<-subs(aus100, data.frame(gridCount100[,c("Var1","val")]))
aus100<-raster::mask(aus100,aus100Mask)
aus50Mask<-subs(aus50, data.frame(gridCount50[,c("Var1","val")]))
aus50<-raster::mask(aus50,aus50Mask)


#extract the gridCell ID for observations
occ$GridCell50ID<-extract(aus50,cbind(occ$Lon_Albers,occ$Lat_Albers),method='simple')#extract the 50km resolution cell value
occ$GridCellID100<-extract(aus100,cbind(occ$Lon_Albers,occ$Lat_Albers),method='simple')#extract the 100km resolution cell value

#############################
#
# #remove obs most likely to be erroneous
#
# #############################
# species<-unique(obsall$Scientific.Name)#names of species
# cleandat<-list()#empty list to feed species dat into
# for(i in 1:length(species)){
#   spdat<-subset(obsall,Scientific.Name==species[i])#get species data
#   sp<-na.omit(unique(data.frame(spdat$GridCellID,as.numeric(1))))#unique grid cells
#   if(nrow(sp)<=3){
#     next
#   }
#   colnames(sp)<-c("ID","loc")
#   splocs<-subs(aus,sp[,c("ID","loc")])# rasterize data
#   #plot(splocs)#plot sp data
#   spPoints<-rasterToPoints(splocs)#turn raster into spatial points dataframe
#   d <- dist1(spPoints) #calculate pairwise distances between points
#   newdata <- as.data.frame(cbind(spPoints, apply(d, 1, function(x) sort(x, decreasing=F)[2])))
#   #Construct new data frame with desired variables
#   colnames(newdata) <- c('x', 'y', 'val', 'distance')
#   nearDist<-mean(newdata$distance)+sd(newdata$distance)+sd(newdata$distance)
#   bad<-subset(newdata,distance>nearDist)#pull out the ones more than 2sd away from mean dist
#   good<-subset(newdata,distance<=nearDist)#get the good ones
#   good$GridCellID<-extract(aus,cbind(good$x,good$y),method='simple')
#   #keep those grid cells with more 4 or more observations
#   bad$GridCellID<-extract(aus,cbind(bad$x,bad$y),method='simple')#
#   spdatbad<-spdat[spdat$GridCellID %in% bad$GridCellID,]
#   keep<-subset(as.data.frame(table(spdatbad$GridCellID)),Freq>=4)
#   if(nrow(keep)>=1){
#     colnames(keep)<-c("GridCellID","val")
#     keep$val<-1
#     #put the data back together and write to list
#     df<-rbind(good[,c("GridCellID","val")],keep)
#   } else{
#     df<-good[,c("GridCellID","val")]
#   }
#
#   df<-spdat[spdat$GridCellID %in% df$GridCellID,]
#   cleandat[[i]]<-df
# }
# obsall<-do.call("rbind",cleandat)
#

######################################

#make plot of species richness and sampling effort 50k

#####################################
effort<-raster::subs(aus50,
                     as.data.frame(table(occ$GridCell50ID)))#raster of sampling effort
cells<-unique(occ[,c("Scientific.Name", "GridCell50ID")])#species in each gridcell
obscountall<-as.data.frame(table(cells$GridCell50ID))#observations in gridcells
richness<-raster::subs(aus50,obscountall)#replace values in raster with number of species

#make plot of species richness and sampling effort 100km
effort100<-raster::subs(aus100,
                        as.data.frame(table(occ$GridCellID100)))#raster of sampling effort
cells100<-na.omit(unique(occ[,c("Scientific.Name", "GridCellID100")]))#species in each gridcell
obscountall100<-as.data.frame(table(cells100$GridCellID100))#observations in gridcells
richness100<-raster::subs(aus100,obscountall100)#replace values in raster with number of species



#########################################

#plot the data to have a look at it

#########################################

par(mfrow=c(1,2),mar=c(2,2,2,0))
plot(effort,
     breaks = c(0,20,100,500,1000,5000,93000),
     col = brewer.pal(7,"BuPu"),
     main="Sampling effort",
     legend=FALSE)
legend(x='topright',
       legend = c("<20","20-100", "100-500", "500-1000","1000-5000",">5000"),
       fill = brewer.pal(7,"BuPu"),
       cex=.5)

#plot the number of species in a gridcell
plot(richness,
     breaks = c(0,10,50,100,200,300,400),
     col = brewer.pal(7,"YlGnBu"),
     main="Species richness",
     legend=FALSE)
legend(x='topright',
       legend = c("<10", "10-50","50-100", "100-200", "200-300",">300"),
       fill = brewer.pal(7,"YlGnBu"),
       cex=.5)



##############################

#use vegan to calculate Extrapolated Species Richness

##############################


#Vexp<-specpool(SpeciesCounts[,2:ncol(SpeciesCounts)],SpeciesCounts$ID, smallsample = TRUE)
#Vexp<-subset(Vexp,n==1)

#rearrange data to use 10km data
occ$lat2<-round(occ$Lat_Albers,5000)
occ$lon2<-round(occ$Lon_Albers,5000)
occ$latlon<-paste(occ$lat2,occ$lon2)#name to group data
SCountLatLon<-as.data.frame(unique(occ$latlon))#10kmgrid name
colnames(SCountLatLon)<-"ID"

species<-unique(occ$Scientific.Name)

for(i in 1:length(species)){
  spdat<-subset(occ,Scientific.Name==species[i])#get species data
  count<-as.data.frame(table(spdat$latlon))#cound in each of the 10km grids
  colnames(count)<-c("ID",species[i])#rename columns
  SCountLatLon<-merge(SCountLatLon,count, by="ID",all.x=TRUE)#merge the species data back with gridcell IDs
  message(i)
}
SCountLatLon[is.na(SCountLatLon)] <- 0 #change NA to 0
lat<-sapply(strsplit(as.vector(SCountLatLon$ID),"\\ "),"[[", 1)#get lats
lon<-sapply(strsplit(as.vector(SCountLatLon$ID),"\\ "),"[[", 2)#get longs
Group<-extract(aus50,cbind(as.numeric(lon),as.numeric(lat)),method='simple')#extract cell ID 50km
Group100<-extract(aus100,cbind(as.numeric(lon),as.numeric(lat)),method='simple')#extract cell ID 100km
########instidence based estimates
Vexp50<-specpool(SCountLatLon[,2:ncol(SCountLatLon)],Group, smallsample = TRUE)#50km
Vexp50$GridCellID<-rownames(Vexp50)

Vexp100<-specpool(SCountLatLon[,2:ncol(SCountLatLon)],Group100, smallsample = TRUE)#100km
Vexp100$GridCellID<-rownames(Vexp100)

###########################

######### Make into plots

###########################

#50km

par(mfrow=c(1,3),mar=c(2,2,2,5)) #Choa estimate at 50km
ER50<-Vexp50[,c("GridCellID","chao")]
ER50$GridCellID<-factor(ER50$GridCellID)
colnames(ER50)<-c("Var1","Freq")
Choa50<-raster::subs(aus50,ER50)#raster
plot(Choa50 ,
     main = "Choa estimate of sp. richness",
     breaks = c(0,50,100,150,200,250,300,400,8000),
     col = c("white",
             'grey87',
             'grey77',
             'grey65',
             'grey45',
             'grey35',
             'grey29',
             'black'),
     legend=FALSE)
legend(x='topright',
       legend = c("<50", "100","150", "200","250","300","400",">400"),
       fill = c("white",
                'grey87',
                'grey77',
                'grey65',
                'grey45',
                'grey35',
                'grey29',
                'black'),
       cex=.5)

Vexp50$freq<-as.numeric(as.character(Vexp50[,"Species"]))#observed species richness
Rich50<-Vexp50[,c("GridCellID","freq")]
Rich50$GridCellID<-factor(Rich50$GridCellID)
colnames(Rich50)<-c("Var1","Freq")
Richness50<-raster::subs(aus50,Rich50)#raster
plot(Richness50,
     main = "Observed species richness", breaks = c(0,50,100,150,200,250,300,400,8000),
     col = c("white",
             'grey87',
             'grey77',
             'grey65',
             'grey45',
             'grey35',
             'grey29',
             'black'),
     legend=FALSE)
legend(x='topright',
       legend = c("<50", "100","150", "200","250","300","400",">400"),
       fill = c("white",
                'grey87',
                'grey77',
                'grey65',
                'grey45',
                'grey35',
                'grey29',
                'black'),
       cex=.5)
#completness
Complete50<-Richness50/Choa50
plot(Complete50,
     main = "Completeness index",
     breaks = c(.2,.4,.6,.8,1),
     col = c("white",
             'grey65',
             'grey29',
             'black'),#brewer.pal(4,"Set1")
     legend=FALSE)
legend(x='topright',
       legend = c("<.4", ".4-.6",".6-.8", ">.8"),
       fill = c("white",
                'grey65',
                'grey29',
                'black'),
       cex=.5)

C50<-Complete50
C50[C50[]>=0]=1
cellStats(C50,stat="sum") #number of cells
C50<-Complete50
C50[C50[]>=.7]=1
C50[C50[]<.7]=NA
plot(C50)
cellStats(C50,stat="sum")#number of cells above .7


#100km
pdf(file = "/Users/daisy/Google Drive/PhD/ThermalStress/figures/ChaoCompletnessIndex20160627.pdf",
     width = 10, height = 3)

par(mfrow=c(1,3),mar=c(2,2,2,3)) #Choa estimate at 100km
ER100<-Vexp100[,c("GridCellID","chao")]
ER100$GridCellID<-factor(ER100$GridCellID)
colnames(ER100)<-c("Var1","Freq")
Choa100<-raster::subs(aus100,ER100)#raster
plot(Choa100 ,
     main = "Chao estimate of species richness",
     breaks = c(0,30,50,70,90,110,130,150,180),
     col = c("white",
             'grey87',
             'grey77',
             'grey65',
             'grey45',
             'grey35',
             'grey29',
             'black'),
     legend=FALSE)
legend(x='topright',
       legend = c("<30", "50","70", "90","110","130","150",">180"),
       fill = c("white",
                'grey87',
                'grey77',
                'grey65',
                'grey45',
                'grey35',
                'grey29',
                'black'),
       cex=.7,
       bty="n")



Vexp100$freq<-as.numeric(as.character(Vexp100[,"Species"]))#observed species richness
Rich100<-Vexp100[,c("GridCellID","freq")]
Rich100$GridCellID<-factor(Rich100$GridCellID)
colnames(Rich100)<-c("Var1","Freq")
Richness100<-raster::subs(aus100,Rich100)#raster
plot(Richness100,
     main = "Observed species richness", breaks = c(0,30,50,70,90,110,130,150,160),
     col = c("white",
             'grey87',
             'grey77',
             'grey65',
             'grey45',
             'grey35',
             'grey29',
             'black'),
     legend=FALSE)
legend(x='topright',
       legend = c("<30", "50","70", "90","110","130","150",">160"),
       fill = c("white",
                'grey87',
                'grey77',
                'grey65',
                'grey45',
                'grey35',
                'grey29',
                'black'),
       cex=.7,
       bty="n")

#completness
Complete100<-Richness100/Choa100

plot(Complete100,
     main = "Completeness index",
     breaks = c(0,.7,1),
     col =  c("white",
              'black'),
     legend=FALSE)
legend(x='topright',
       legend = c("< 0.7","> 0.7"),
       fill = c("white",
                'black'),
       cex=.7,bty = "n")


dev.off()



#how many grid cells less than >.
C100<-Complete100
C100[C100[]>=.7]=1
C100[C100[]<.7]=NA
plot(C100)

writeRaster(trim(C100),
            paste0("/Users/daisy/Google Drive/PhD/ThermalStress/data/ChoaIndex100km",as.Date(Sys.time()),".asc"),
            overwrite=TRUE)


##################

#Remove grid cells where there are less than 10 domed species or 10 cupped species
#read in trait database, this includes nest type information
traits<-read.csv('/Users/daisy/Google Drive/PhD/ThermalStress/tables/Table_S1_20161128.csv')


domeSP<-subset(traits,Nest.type=="Dome")
domeOcc<-occ[occ$Scientific.Name %in% domeSP$Species,]
domeOcc<-unique(domeOcc[,c("Scientific.Name","GridCellID100")])
domecount<-as.data.frame(table(domeOcc$GridCellID100))
domeRich<-raster::subs(aus100,domecount)
domeRich[domeRich[]<10]=NA
domeRich<-domeRich/domeRich

cupSP<-subset(traits,Nest.type=="Cup")
cupOcc<-occ[occ$Scientific.Name %in% cupSP$Species,]
cupOcc<-unique(cupOcc[,c("Scientific.Name","GridCellID100")])
cupcount<-as.data.frame(table(cupOcc$GridCellID100))
cupRich<-raster::subs(aus100,cupcount)
cupRich[cupRich[]<10]=NA
cupRich<-cupRich/cupRich

endComplete<-C100+domeRich+cupRich
endComplete[endComplete[]<3]=NA
endComplete<-endComplete/endComplete
writeRaster(trim(endComplete),
            paste0("/Users/daisy/Google Drive/PhD/ThermalStress/data/ChoaIndex100km10SpeciesMinimum",as.Date(Sys.time()),".asc"),
            overwrite=TRUE)

