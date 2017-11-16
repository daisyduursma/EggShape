

#clean and combine breeding observation data from all museums, ATLAS, Nest Record Scheme, ebird, ala, gbif. 
#This is limited to Passarines
#fix lat. values that are not negative
#check observations fall within Australia based on AWAP data (same data that is planned to use on thermal stress study)



rm(list = ls())

library(raster)
#library(car) 
#library(stringr)
library(Hmisc)
#library(lubridate)
library(gtools)
library(data.table)

indir<-'~/Google Drive/PhD/Data/Observaitons/Raw/egg'

# gnt<-subset(read.csv('~/Google Drive/PhD/birdTraits/PublishedData/Australian_Bird_Data_Version_1.csv'),
#             X31_Non.breeding_populations_of_core_taxa_4==0 &
#               X32_Extinct_4==0 &
#               X34_Vagrant_4==0 &
#               X17_Species_2==1)                                     
# 
# gnt<-gnt[c("Species","X3_Taxon_common_name_2")]

traits<-read.csv('/Users/daisy/Google Drive/PhD/ThermalStress/data/ThermalStressTraits2016-05-09.csv')
traits<-subset(traits,USEThermalStressStudy != 'omit')

species<-traits$Species
JetzSpecies<-traits$Jetz.name


# 
# #mask of Australia
# aus<-raster("~/Google Drive/PhD/Data/Spatial/Climate/EMASTBiovars/bio_1.asc")
# aus_XY <- as.data.frame(rasterToPoints(aus))#find lat and long of 1km data
# aus<-aggregate(aus,fact=50,na.rm=TRUE)#make raster 1 degree
# aus[]<-1:5810 #Give each grid cell unique value
# 
# #keep only grid cells that are at least 1/3 land
# aus_XY$GridCellID<-extract(aus,cbind(aus_XY$x,aus_XY$y),method='simple')#extract the 1km resolution cell value
# gridCount<-as.data.frame(table(aus_XY$GridCellID))#find how many grid cells fall in the 1 degree data
# gridCount<-subset(gridCount,Freq>=833)#only keep grid cells that are actually at least 33% land
# gridCount$val<-1 #set the value to 1
# 
# #make a mask where land is at least 33% of the 1 degree grid cells
# ausMask<-subs(aus, data.frame(gridCount[,c("Var1","val")])) #make a 1 degree mask where gridcells are at least 33% land
# plot(ausMask)
# rm(aus_XY)


###########################################################################
######NT data
NT<-read.csv(paste0(indir,'/NT_Museum/MAGNT_bird eggs_Jan2015.csv'))
NT<-subset(NT, !is.na(Lat_degrees) & !is.na(Long_degrees))
NT$Lat_minutes <- replace(NT$Lat_minutes, is.na(NT$Lat_minutes),0)
NT$Long_minutes <- replace(NT$Long_minutes, is.na(NT$Long_minutes),0)
NT$lat<-round(NT$Lat_degrees + (NT$Lat_minutes/60),2)
NT$long<-round(NT$Long_degrees + (NT$Long_minutes/60),2)
NT$Sci.Name<-paste(NT$Genus,NT$Species,sep=' ')
eggdat<-NT[c('Sci.Name','lat','long')]
colnames(eggdat)<-c( 'Scientific.Name','lat','lon')

###########################################################################
######South Australia Museaum data
SA<-read.csv(paste0(indir,'/SA_Museam/SA_egg_records.csv'))
SA$lat_min <- replace(SA$lat_min, is.na(SA$lat_min),0)
SA$long_min <- replace(SA$long_min, is.na(SA$long_min),0)
SA<-subset(SA, !is.na(lat_hour) & !is.na(lat_hour))
SA$lat<-round(SA$lat_hour + (SA$lat_min/60),2)
SA$long<-round(SA$long_hour + (SA$long_min/60),2)
SA$Sci.Name<-SA$Taxon
SA<-SA[c('Sci.Name','lat','long')]
colnames(SA)<-c( 'Scientific.Name','lat','lon')

eggdat<-smartbind(eggdat,SA)
rm(SA)

###########################################################################
######Queensland
QU<-read.csv(paste0(indir,'/QL_Museum/egg_nest_sepLATLONG.csv'))
# QU<-subset(QU,Field.Coll.Specimen.Category=="Egg(s)" |
#              Field.Coll.Specimen.Category=="Egg(s), Nest" |
#              Field.Coll.Specimen.Category== "Egg(s), Skeletal parts Nest" |
#              Field.Coll.Specimen.Category== "Spirit, Egg(s)")
QU$lat_min <- replace(QU$lat_min, is.na(QU$lat_min),0)
QU$long_min <- replace(QU$long_min, is.na(QU$long_min),0)
QU$lat_sec <- replace(QU$lat_sec, is.na(QU$lat_sec),0)
QU$Long_sec <- replace(QU$Long_sec, is.na(QU$Long_sec),0)
QU$lat<-round(QU$lat_deg + (QU$lat_min/60) + (QU$lat_sec/3600) *-1,2)
QU$long<-round(QU$long_deg + (QU$long_min/60)+ (QU$Long_sec/3600) *-1,2)
QU$Sci.Name<-paste(QU$Genus, QU$Species, sep=" ")
QU<-QU[c('Sci.Name','lat','long')]
colnames(QU)<-c( 'Scientific.Name','lat','lon')

eggdat<-smartbind(eggdat,QU)

rm(QU)
###########################################################################
######Australia Museaum data
ausmus<-read.csv(paste0(indir,'/Australia Museum/AUSMUS_egg_records.csv'))
ausmus<-subset(ausmus, !is.na(start_latitude) & !is.na(start_longitude))
ausmus$Sci.Name<-paste(ausmus$genus,ausmus$species,sep=' ')
ausmus<-ausmus[c('Sci.Name','start_latitude','start_longitude')]
colnames(ausmus)<-c( 'Scientific.Name','lat','lon')
ausmus$lat<-round(ausmus$lat,2)
ausmus$lon<-round(ausmus$lon,2)
eggdat<-smartbind(eggdat,ausmus)
rm(ausmus)
###########################################################################
######Victoria

Vic1<-read.csv(paste0(indir,'/Vic_Museum/OZCAM/data.csv'))
Vic1<-Vic1[,c('Scientific.Name','Latitude...processed','Longitude...processed')]
colnames(Vic1)<-c('Scientific.Name','lat', 'lon')
Vic2<-read.csv(paste0(indir,'/Vic_Museum/MVEgg Collection1.csv'))[c('Genus...Version.1.2.elements..1..', 'Species...Version.1.2.elements..2..','Latitude...Version.1.2.elements..3..','Longitude...Version.1.2.elements..3..')]
Vic3<-read.csv(paste0(indir,'/Vic_Museum/MVEgg Collection2.csv'))[c('Genus...Version.1.2.elements..1..', 'Species...Version.1.2.elements..2..','Latitude...Version.1.2.elements..3..','Longitude...Version.1.2.elements..3..')]
Vic4<-read.csv(paste0(indir,'/Vic_Museum/MVEgg Collection3.csv'))[c('Genus...Version.1.2.elements..1..', 'Species...Version.1.2.elements..2..','Latitude...Version.1.2.elements..3..','Longitude...Version.1.2.elements..3..')]
Vic5<-read.csv(paste0(indir,'/Vic_Museum/MVEgg Collection4.csv'))[c('Genus...Version.1.2.elements..1..', 'Species...Version.1.2.elements..2..','Latitude...Version.1.2.elements..3..','Longitude...Version.1.2.elements..3..')]
NVic<-rbind(Vic2,Vic3,Vic4,Vic5)
NVic$Species<-paste(NVic$Genus...Version.1.2.elements..1..,NVic$'Species...Version.1.2.elements..2..')
NVic<-subset(NVic, !is.na('Latitude...Version.1.2.elements..3..') 
             & !is.na('Longitude...Version.1.2.elements..3..')
             & !is.na('Scientific.Name'))
NVic<-NVic[c('Species','Latitude...Version.1.2.elements..3..','Longitude...Version.1.2.elements..3..')]
colnames(NVic)<-c( 'Scientific.Name','lat','lon')
NVic<-rbind(Vic1,NVic)
NVic$lon<-round(NVic$lon,2)
NVic$lat<-round(NVic$lat,2)
eggdat<-na.omit(smartbind(eggdat,NVic))
rm(Vic1,Vic2,Vic3,Vic4,Vic5,NVic)
###########################################################################
######ANWC
ANWC<-read.csv(paste0(indir,'/ANWC_OZCAM/OZCAM_eggs.csv'))
ANWC<-subset(ANWC,Country...parsed=='Australia')
ANWC<-subset(ANWC, !is.na(Latitude...processed) 
             & !is.na(Longitude...processed))
ANWC<-ANWC[c('Scientific.Name','Latitude...processed','Longitude...processed')]
colnames(ANWC)<-c( 'Scientific.Name','lat','lon')
ANWC$lon<-round(ANWC$lon,2)
ANWC$lat<-round(ANWC$lat,2)
eggdat<-smartbind(eggdat,ANWC)
rm(ANWC)

###########################################################################
###### Western Australia
WA<-read.csv(paste0(indir,'/WA_Museum/WA_museum.csv'))
WA$precision<-ifelse(!is.na(WA$Coordinate.Uncertainty.in.Metres...parsed),
                     WA$Coordinate.Uncertainty.in.Metres...parsed/1000,"NA")
#get the ID for egg records
#rego<-read.csv(paste0(indir,'/WA_Museum/eggs-regno.csv'))[,1]
WAeggs<-WA
WAeggs<-subset(WAeggs, Country...parsed=='Australia')
WAeggs<-subset(WAeggs, !is.na(Latitude...processed) 
               & !is.na(Longitude...processed))
WAeggs<-WAeggs[c('Scientific.Name','Latitude...processed','Longitude...processed')]
colnames(WAeggs)<-c( 'Scientific.Name','lat','lon')
WAeggs$lon<-round(WAeggs$lon,2)
WAeggs$lat<-round(WAeggs$lat,2)
eggdat<-smartbind(eggdat,WAeggs)

rm(WAeggs,WA)

###########################################################################
##### Nest Record Scheme
NRS<-read.csv(paste0(indir,'/BLA_NRS/NRSExtract.csv'))
NRS<-with(NRS,NRS[!is.na(Lat) & !is.na(Lon),])
NRS<-NRS[c('Scientific_name','Lat','Lon')]
colnames(NRS)<-c( 'Scientific.Name','lat','lon')
NRS$lon<-round(NRS$lon,2)
NRS$lat<-round(NRS$lat,2)
eggdat<-smartbind(eggdat,NRS)
rm(NRS)


###########################################################################
######### Queen Victoria Museum 
QV<-read.csv(paste0(indir,'/QV_Museum/QVM.csv'))
QV<-QV[,c('Scientific.Name','Latitude...processed','Longitude...processed')]
colnames(QV)<-c( 'Scientific.Name','lat','lon')
QVegg<-read.csv(paste0(indir,'/QV_Museum/QVM Egg Data.csv'))
QVegg$Scientific.Name<-paste(QVegg$Genus, QVegg$Species,sep = ' ')
QVegg<-QVegg[,c('Scientific.Name','lat','long')]
colnames(QVegg)<-c( 'Scientific.Name','lat','lon')
QV<-rbind(QV,QVegg)
QV<-subset(QV, !is.na(lat) & !is.na(lon))
QV<-QV[c('Scientific.Name','lat','lon')]
colnames(QV)<-c( 'Scientific.Name','lat','lon')
QV$lon<-round(QV$lon,2)
QV$lat<-round(QV$lat,2)
eggdat<-smartbind(eggdat,QV)
rm(QV,QVegg)

##############################################################################
#TMAG
tmag<-read.csv(paste0(indir,'/TMAG_Museum/TMAGeggsnests.csv'))
tmag<-subset(tmag, !is.na(Lat) & !is.na(Long))
tmag$date<-paste0(tmag$year,'-',tmag$month,'-',tmag$day)
tmag<-tmag[c('Species.Name','Lat','Long')]
colnames(tmag)<-c( 'Scientific.Name','lat','lon')
tmag$lon<-round(tmag$lon,2)
tmag$lat<-round(tmag$lat,2)
eggdat<-smartbind(eggdat,tmag)
rm(tmag)
eggdat<- eggdat[!duplicated(eggdat),]#remove duplicates

eggdat$Breeding<-1
##############################################################################
#ATLAS
ATLAS<-fread('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/observation/ATLAS/AtlasSurveysObs.csv',data.table=FALSE)
#add Scientific names
birdNames<-read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/BLA_Working_List_v1.1.csv')
ATLAS<-merge(ATLAS,birdNames,by.x = 'AtlasNo',by.y='SpNo',all.x=TRUE)
ATLAS<-subset(ATLAS,!is.na(Taxon.scientific.name))
ATLAS<-ATLAS[,c('Taxon.scientific.name','Lat','Lon','Breeding')]
colnames(ATLAS)<-c( 'Scientific.Name','lat','lon','Breeding')
ATLAS$Breeding<-ifelse(ATLAS$Breeding!=1,NA,1)
#replace 0 with NA
ATLAS$Breeding[is.na(ATLAS$Breeding)] <- 0
ATLAS<- na.omit(ATLAS[!duplicated(ATLAS),])
eggdat<-smartbind(eggdat,ATLAS)

eggdat$lat<-round(eggdat$lat, digits = 2)
eggdat$lon<-round(eggdat$lon, digits = 2)
eggdat<- eggdat[!duplicated( eggdat),]

rm(ATLAS)


##################################################################
#########ABBBS - keep only one unique lat, long, day


abbbsA<-read.csv(paste0(indir,"/ABBBS/pullusBandingRecordsPartA_dd.csv"))
abbbsB<-read.csv(paste0(indir,"/ABBBS/pullusBandingRecordsPartB_dd.csv"))
abbbs<-rbind(abbbsA,abbbsB)
abbbs<-abbbs[!duplicated(abbbs[,c("SCIENTIFIC_NAME","LAT","LON")]),]
abbbs$lat<-round(abbbs$LAT,2)
abbbs$lon<-round(abbbs$LON,2)
abbbs$Scientific.Name<-abbbs$SCIENTIFIC_NAME 
abbbs$Breeding<-1
abbbs<-abbbs[,c('Scientific.Name','lat','lon','Breeding')]
eggdat<-smartbind(eggdat,abbbs)

rm(abbbs,abbbsB,abbbsA)

#######################
#tidy up before continuing
#keep only gridcell ID to reduce number
latfix<-subset(eggdat, lat > 0)
latfix$lat<-latfix$lat*-1
latgood<-subset(eggdat, lat <= 0)
eggdat<-rbind(latfix, latgood)
rm(latfix,latgood)
eggdat<-eggdat[!duplicated(eggdat),]

################
#ALA

aves2 <- read.delim("~/Google Drive/PhD/Data/Observaitons/Raw/observation/ALA/aves2.csv")
aves2<-aves2[,c("taxon_name","latitude","longitude")]
aves2$latitude<-round(aves2$latitude, digits = 2)
aves2$longitude<-round(aves2$longitude, digits = 2)
aves2<- na.omit(aves2[!duplicated( aves2),])
colnames(aves2)<-c( 'Scientific.Name','lat','lon')
aves2$Breeding<-0


eggdat<-smartbind(eggdat,aves2)
rm(aves2)
eggdat<- eggdat[!duplicated(eggdat),]

##################################################################
#########eBird
ebird<-read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/observation/eBird/ebd_AU_prv_relFeb-2015/ebd_AU_prv_relFeb-2015.txt',
                header=TRUE,sep="\t",strip.white=TRUE)
ebird$lat<-round(ebird$LATITUDE,2)
ebird$lon<-round(ebird$LONGITUDE,2)
ebird$Scientific.Name<-ebird$SCIENTIFIC.NAME
ebird$Breeding<-with(ebird,ifelse (BREEDING.BIRD.ATLAS.CODE =='NY'|BREEDING.BIRD.ATLAS.CODE =='NE'|
                                 BREEDING.BIRD.ATLAS.CODE =='ON'|BREEDING.BIRD.ATLAS.CODE =='PE'|
                                 BREEDING.BIRD.ATLAS.CODE =='FL'|BREEDING.BIRD.ATLAS.CODE =='CS' |
                                 BREEDING.BIRD.ATLAS.CODE =='FY'|BREEDING.BIRD.ATLAS.CODE =='N',1,0))
ebird<-ebird[c('Scientific.Name','lat','lon','Breeding')]
ebird$lat<-round(ebird$lat,2)
ebird$lon<-round(ebird$lon,2)
ebird<- na.omit(ebird[!duplicated( ebird),])
eggdat<-smartbind(eggdat,ebird)
eggdat<- eggdat[!duplicated(eggdat),]
rm(ebird)

#######################
gbif<-fread("~/Google Drive/PhD/Data/Observaitons/Raw/observation/Gbif/0015069-160118175350007.csv",data.table=FALSE)
gbif$lat<-round(gbif$decimallatitude,2)
gbif$lon<-round(gbif$decimallongitude,2)
gbif$Scientific.Name<-gbif$species
gbif<-gbif[c('Scientific.Name','lat','lon')]
gbif<- na.omit(gbif[!duplicated( gbif),])
gbif$Breeding<-0
eggdat<-smartbind(eggdat,gbif)
eggdat<- eggdat[!duplicated(eggdat),]
eggdat$lat<-round(eggdat$lat,2)
eggdat$lon<-round(eggdat$lon,2)
eggdat<- eggdat[!duplicated(eggdat),]

rm(gbif)

#keep only gridcell ID to reduce number
latfix<-subset(eggdat, lat > 0)
latfix$lat<-latfix$lat*-1
latgood<-subset(eggdat, lat <= 0)
eggdat<-rbind(latfix, latgood)
rm(latfix,latgood)


#################################################################
# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#finish cleaning egdat
eggdat$Scientific.Name<- capitalize(trim(gsub('  ', ' ', 
                                              eggdat$Scientific.Name)))#check for extra spaces and capatlize
#remove duplicates
eggdat<-eggdat[!duplicated(eggdat),]


#fix up names
syn<-subset(read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/namesResolved20160224.csv'), use ==1)
bad<-subset(read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/namesResolved20160224.csv'), use ==0)
bad$scn1<- capitalize(trim(gsub('  ', ' ', 
                                bad$scn1)))#check for extra spaces and capatlize
syn$scn1<- capitalize(trim(gsub('  ', ' ', 
                                syn$scn1)))#check for extra spaces and capatlize

syn$garnett_name<- capitalize(trim(gsub('  ', ' ', 
                                syn$garnett_name)))

eggdat$Scientific.Name<-gsub("[ ]\\(.+\\)","",eggdat$Scientific.Name)
obsSmall<-eggdat[!duplicated(eggdat),]
rm(eggdat)

obsSmall<-obsSmall[obsSmall$Scientific.Name %nin% bad$scn1,]
fix<-obsSmall[obsSmall$Scientific.Name %in% syn$scn1,]
nofix<-obsSmall[obsSmall$Scientific.Name %nin% syn$scn1,]
fix<-merge(fix,syn,by.x = 'Scientific.Name',by.y='scn1',all.x=TRUE)
fix<-fix[!duplicated(fix),]
fix$Scientific.Name <-fix$garnett_name
fix<-fix[,c("Scientific.Name","lat","lon","Breeding")]
obsSmall<-smartbind(fix,nofix)
obsSmall<-na.omit(obsSmall[!duplicated(obsSmall),])

#make sure data is within continental Australia
aus<-raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/BOM/maxann.txt",
           crs = '+proj=longlat +datum=WGS84')

finalobsSmall<-obsSmall[obsSmall$Scientific.Name %in% species,]
finalobsSmall2<-obsSmall[obsSmall$Scientific.Name %in% JetzSpecies,]
finalobsSmall<-as.data.frame(rbind(finalobsSmall,finalobsSmall2))
finalobsSmall<-na.omit(finalobsSmall[!duplicated(finalobsSmall),])


#make sure data is within continental Australia
aus<-raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/BOM/maxann.txt",
crs = '+proj=longlat +datum=WGS84')

xy<-cbind(finalobsSmall$lon,finalobsSmall$lat)
finalobsSmall$outsideAustralia<-is.na(extract(aus,xy))
finalobsSmall<-subset(finalobsSmall,outsideAustralia==FALSE)[,c( 'Scientific.Name','lat','lon','Breeding')]

write.csv(finalobsSmall, paste0('~/Google Drive/PhD/Data/Observaitons/Cleaned/observations/PASSARINE_ALAGBIFABBBSATLASGridCellIDObservations',
                              as.Date(Sys.time()),'.csv'),row.names=FALSE)

