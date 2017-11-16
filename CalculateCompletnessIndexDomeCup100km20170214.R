#calculate completeness index 100km gridcells
#and remove cells with not enough data 
#proform on open and closed nest species

rm(list = ls())
library(raster)
library(Hmisc)
library(data.table)
library(rgeos)
library(SpatialTools)
library(RColorBrewer)
library(gtools)
library(vegan)

#data for equal area
Albers<-"+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
km1<-raster(nrows=4900, ncols=5600, 
            xmn=-2904860, xmx=2695140,ymn=-5334205, ymx=-434205,
            crs =paste(Albers))

#read in observations
occ <- fread("~/Google Drive/PhD/ThermalStress/data/PassarineObservationsFinal2016-06-21.csv",
                        data.table=FALSE)
#change XY to albers equal area
Locs_sp <- SpatialPoints(cbind(occ$lon, occ$lat), proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
Locs_sp <- spTransform(Locs_sp, CRS=CRS(paste(Albers)))
Alberslocs<-as.data.frame(Locs_sp)
colnames(Alberslocs) <- c("Lon_Albers","Lat_Albers")
occ<-cbind(occ,Alberslocs)

species <- unique(occ$Scientific.Name)

#read in trait database, this includes nest type information
traits<-read.csv('/Users/daisy/Google Drive/PhD/ThermalStress/tables/Table_S1_20161128.csv')

aus<-raster::raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/BOM/maxann.txt",
                    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")#mask of Australia
aus<-aus/aus #make it into a maks
ausAlbers<-projectRaster(aus,km1,method='ngb')#project to albers
aus_XY <- as.data.frame(rasterToPoints(ausAlbers))#find lat and long of 1km data
# aus100<-aggregate(ausAlbers,110,'mean')#100kmX100km resolution
# aus100[]<-1:2295 #Give each grid cell unique value
# 
# aus_XY$GridCell100ID<-raster::extract(aus100,cbind(aus_XY$x,aus_XY$y),method='simple')#extract the 1km resolution cell value
# gridCount100<-as.data.frame(table(aus_XY$GridCell100ID))
# gridCount100<-subset(gridCount100,Freq>=4033)#only keep grid cells that are actually at least 33% land
# gridCount100$val<-1 #set the value to 1
# aus100Mask<-subs(aus100, data.frame(gridCount100[,c("Var1","val")]))
# aus100<-raster::mask(aus100,aus100Mask)

aus100<-aggregate(ausAlbers,100,'mean')#100kmX100km resolution
aus100[]<-1:2744 #Give each grid cell unique value

aus_XY$GridCell100ID<-raster::extract(aus100,cbind(aus_XY$x,aus_XY$y),method='simple')#extract the 1km resolution cell value
gridCount100<-as.data.frame(table(aus_XY$GridCell100ID))
gridCount100<-subset(gridCount100,Freq>=3333)#only keep grid cells that are actually at least 33% land
gridCount100$val<-1 #set the value to 1
aus100Mask<-subs(aus100, data.frame(gridCount100[,c("Var1","val")]))
aus100<-raster::mask(aus100,aus100Mask)


#extract the gridCell ID for observations
occ$GridCellID100<-extract(aus100,cbind(occ$Lon_Albers,occ$Lat_Albers),method='simple')#extract the 100km resolution cell value
occ<-subset(occ,!is.na(GridCellID100))

#make raster of sampling effort
effort<-raster::subs(aus100,
                        as.data.frame(table(occ$GridCellID100)))#raster of sampling effort

effort[effort<5000]<-0
effort[effort>=5000]<-1

#round to base number (0.5)
mround <- function(x,base){ 
      base*round(x/base) 
} 

# Simple function for placing labels on a figure.
plotlabel <- function(txt, where, inset=0.2, font=2, inset.x=inset, inset.y=inset,...){
  u <- par()$usr
  if(grepl("left",where))x <- u[1] + inset.x*(u[2]-u[1])
  if(grepl("right",where))x <- u[2] - inset.x*(u[2]-u[1])
  if(grepl("bottom",where))y <- u[3] + inset.y*(u[4]-u[3])
  if(grepl("top",where))y <- u[4] - inset.y*(u[4]-u[3])
  
  text(x,y,txt,font=font,...)
}


ChoaNestType<-function (occurances,ausMASK){
  #rearrange data to use 10km data
  occurances$lat2<-mround(occurances$Lat_Albers,5000)
  occurances$lon2<-mround(occurances$Lon_Albers,5000)
  occurances$latlon<-paste(occurances$lat2,occurances$lon2)#name to group data
  SCountLatLon<-as.data.frame(unique(occurances$latlon))#5kmgrid name
  colnames(SCountLatLon)<-"ID"
  species<-unique(occurances$Scientific.Name)
  
  for(i in 1:length(species)){
    spdat<-subset(occurances,Scientific.Name==species[i])#get species data
    count<-as.data.frame(table(spdat$latlon))#cound in each of the 10km grids
    colnames(count)<-c("ID",species[i])#rename columns
    SCountLatLon<-merge(SCountLatLon,count, by="ID",all.x=TRUE)#merge the species data back with gridcell IDs
    message(i)
  }
  SCountLatLon[is.na(SCountLatLon)] <- 0 #change NA to 0
  lat<-sapply(strsplit(as.vector(SCountLatLon$ID),"\\ "),"[[", 1)#get lats
  lon<-sapply(strsplit(as.vector(SCountLatLon$ID),"\\ "),"[[", 2)#get longs
  Group100<-extract(ausMASK,cbind(as.numeric(lon),as.numeric(lat)),method='simple')#extract cell ID 100km
  
  ########instidence based estimates
  Vexp100<-specpool(SCountLatLon[,2:ncol(SCountLatLon)],Group100, smallsample = TRUE)#100km
  Vexp100$GridCellID<-rownames(Vexp100)
  return(Vexp100)
}  



#domed
domeSP<-subset(traits,Nest.type=="Dome")
domeOcc<-occ[occ$Scientific.Name %in% domeSP$Species,]
DomedChoa<-ChoaNestType(domeOcc,aus100)#run the choa function
#Domed Rasters
  #Estimted Richness
    ERdome<-DomedChoa[,c("GridCellID","chao")]
    ERdome$GridCellID<-factor(ERdome$GridCellID)
    colnames(ERdome)<-c("Var1","Freq")
    ChoaDomeR<-raster::subs(aus100,ERdome) # raster Choa estimated Richness
  #Obserced Richness
    DomedChoa$freq<-as.numeric(as.character(DomedChoa[,"Species"]))
    RichDome<-DomedChoa[,c("GridCellID","freq")]
    RichDome$GridCellID<-factor(RichDome$GridCellID)
    colnames(RichDome)<-c("Var1","Freq")
    RichDomeR<-raster::subs(aus100,RichDome)#raster
    RichDomeR[RichDomeR<10]<-NA#mask out gridcells with less than 10 species observed
  #completeness Index
    CmpltDomed<-RichDomeR/ChoaDomeR
    CmpltDomed[CmpltDomed<.7]<-NA
    CmpltDomed[CmpltDomed>=.7]<-1
  #add in gridcells whereRichenss is above 10 and effort is greater than 5000  
    dome10<-RichDomeR
    dome10[dome10>=10]<-1
    dome10<-dome10+effort
    dome10[dome10<2]<-NA
    dome10[dome10==2]<-1
    dome10[is.na(dome10)]<-0
    CmpltDomed[is.na(CmpltDomed)]<-0
    CmpltDomed<-CmpltDomed+dome10
    CmpltDomed[CmpltDomed==2]<-1
    CmpltDomed[CmpltDomed==0]<-NA
    
#cupped
cupSP<-subset(traits,Nest.type=="Cup")
cupOcc<-occ[occ$Scientific.Name %in% cupSP$Species,]
CupChoa<-ChoaNestType(cupOcc,aus100)#run the choa function
#cupped Rasters
  #Estimted Richness
    ERcup<-CupChoa[,c("GridCellID","chao")]
    ERcup$GridCellID<-factor(ERcup$GridCellID)
    colnames(ERcup)<-c("Var1","Freq")
    ChoaCupR<-raster::subs(aus100,ERcup) # raster Choa estimated Richness
  #Obserced Richness
    CupChoa$freq<-as.numeric(as.character(CupChoa[,"Species"]))
    RichCup<-CupChoa[,c("GridCellID","freq")]
    RichCup$GridCellID<-factor(RichCup$GridCellID)
    colnames(RichCup)<-c("Var1","Freq")
    RichCupR<-raster::subs(aus100,RichCup)#raster
    RichCupR[RichCupR<10]<-NA#mask out gridcells with less than 10 species observed
  #completeness Index
    CmpltCup<-RichCupR/ChoaCupR
    CmpltCup[CmpltCup<.7]<-NA
    CmpltCup[CmpltCup>=.7]<-1
  #add in gridcells whereRichenss is above 10 and effort is greater than 5000  
    cup10<-RichCupR
    cup10[cup10>=10]<-1
    cup10<-cup10+effort
    cup10[cup10<2]<-NA
    cup10[cup10==2]<-1
    cup10[is.na(cup10)]<-0
    CmpltCup[is.na(CmpltCup)]<-0
    CmpltCup<-CmpltCup+cup10
    CmpltCup[CmpltCup==2]<-1
    CmpltCup[CmpltCup==0]<-NA
    

BWcolD = c("white", 'grey65', 'grey35', 'black')
RichnBrksD = c(0,9.99,20,40,60)
LegD<-c("<10", "20","40","<60")
BWcolC = c("white",'grey65','grey45','grey35','grey20','black')
RichnBrksC = c(0,10,25,50,75,100,200)
LegC<-c("<10", "25","50","75","100",">100")

{
pdf(file = paste0("/Users/daisy/Google Drive/PhD/ThermalStress/figures/ChaoCompletnessIndex",
                  as.Date(Sys.time()),".pdf"),
                  width = 7, height = 4)

par(mfrow=c(2,3), xpd=NA, oma=c(0,1,1,1), mar=c(0,1,1,3), cex=.7, xaxs="i",
    yaxs = "i")
#Choa estimate at 100km

#Domed observed Richenss
plot(RichDomeR,main = "Observed richness",font.main = 1,cex.main=1,axes=FALSE, box=FALSE, breaks = RichnBrksD,col = BWcolD, legend=FALSE)
legend(x='topleft',legend = LegD,fill = BWcolD,bty="n",cex=.8,inset=c(-.02,0))
mtext(expression(paste("(", bold("A"), ") Dome")),line=0,adj=0,cex=.7,outer=TRUE)

#Domed Choa estimated Richenss
plot(ChoaDomeR , main = "Estimated richness",font.main = 1,cex.main=1,axes=FALSE, box=FALSE,breaks = RichnBrksD,col = BWcolD,legend=FALSE)
legend(x='topleft',legend = LegD,fill = BWcolD,bty="n",cex=.8,inset=c(-.02,0))
#Domed completness
plot(CmpltDomed,main = "Completeness index",font.main = 1,cex.main=1,axes=FALSE, box=FALSE,breaks = c(0,.699,1),col =  c("white",'black'),legend=FALSE)
legend(x='topleft',legend = c("< 0.7",">= 0.7"),fill = c("white",'black'),cex=.8,bty = "n",inset=c(-.02,0))


#Cup observed Richenss
plot(RichCupR, axes=FALSE, box=FALSE,breaks = RichnBrksC,col = BWcolC, legend=FALSE)
legend(x='topleft',legend = LegC,fill = BWcolC,cex=.8,bty="n",inset=c(-.02,-.03))
mtext(expression(paste("(", bold("B"), ") Cup")),line=-15,adj=0,cex=.7,outer=TRUE)

#Cup Choa estimated Richenss
plot(ChoaCupR ,axes=FALSE, box=FALSE,breaks = RichnBrksC,col = BWcolC,legend=FALSE)
legend(x='topleft',legend = LegC,fill = BWcolC,cex=.8,bty="n",inset=c(-.02,-.03))
#Cup completness
plot(CmpltCup,axes=FALSE, box=FALSE,breaks = c(0,.699,1),col =  c("white",'black'),legend=FALSE)
legend(x='topleft',legend = c("< 0.7",">= 0.7"),fill = c("white",'black'),cex=.8,bty = "n",inset=c(-.02,-.03))
dev.off()

}

######################################
#write out gridcells that will be used
######################################

CmpltCells<-CmpltCup*CmpltDomed
writeRaster(CmpltCells,
            paste0("/Users/daisy/Google Drive/PhD/ThermalStress/data/CmpltDomeCupGridCells100km",as.Date(Sys.time()),".asc"),
            overwrite=TRUE)



