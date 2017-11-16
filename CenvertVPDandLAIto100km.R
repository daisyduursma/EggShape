#convert VPD and LAI index to Cyclindrical Equal Area 100km resolution

rm(list = ls())
library(raster)
library(rgdal)
library(plantecophys)
library(stringr)
library(RCurl)

#Albers equal area projection
Albers<-"+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#make 100km equal area raster Albers equal area:
km100 <- raster(nrows=49, ncols=56, 
                xmn=-2904860, xmx=2695140,ymn=-5334205, ymx=-434205,
                crs =paste(Albers))
km1<-raster(nrows=4900, ncols=5600, 
            xmn=-2904860, xmx=2695140,ymn=-5334205, ymx=-434205,
            crs =paste(Albers))

#########################
#LAI
########################

samples = 4790
lines   = 3726
file <-"/Users/daisy/Google Drive/PhD/Data/Spatial/LAI/modis_climatology_avg_over_yr_daisy.bin"
dat <- readBin(file,what="double", size=8, n = lines * samples)
LAI <- raster(nrow=lines,ncol=samples)
ext<-extent(c(110.001174,154.998826,-44.998826,-10.001174))
extent(LAI)<-ext
LAI[] <- dat
LAI[LAI<=-999] <- NA

# LAIp <- rasterToPoints(LAI,spatial=TRUE)
# LAIpTF <- spTransform(LAI, CRS=CRS(paste(Albers)))


LAIAlbers1<-projectRaster(LAI,km1,method='bilinear')
LAI100<-aggregate(LAIAlbers1,100,'mean')

writeRaster(LAI100,paste0("/Users/daisy/Google Drive/PhD/ThermalStress/data/spatial/LAI100km",as.Date(Sys.time()),".asc"))



###########################
#VPD
############################


VPDmaker <- function(vapour, maxTemp){
  v<-as.data.frame(rasterToPoints(vapour))
  v$tmp2<-v$tmp2*100 #cnvert to hPa
  Tx<-as.data.frame(rasterToPoints(maxTemp))
  Tx$Vsat<-esat(Tx$tmp) #calculate saturated vapour pressure
  Tx$VPD<-(Tx$Vsat-v$tmp2)/1000 #convert to kPa
  VPD<-rasterFromXYZ(cbind(Tx$x,Tx$y,Tx$VPD))
  return(VPD)
}


#dates
dates<-seq(as.Date("1950/01/01"), as.Date("2015/12/31"), "day")
##download daily data and convert to VPD, write each daily data our because this has to be done ~25000 times
# for(i in 1:length(dates)){
#   D<-strftime(dates[i],format = "%d")
#   D<-str_pad(D,2, pad="0")
#   M<-strftime(dates[i],format = "%m")
#   M<-str_pad(M,2, pad="0")
#   Y<-strftime(dates[i],format = "%Y")
#   
#   temp <- tempfile()#prep temperture data
#   t <- paste0("http://www.bom.gov.au/web03/ncc/www/awap/temperature/maxave/daily/grid/0.05/history/nat/",Y,M,D,Y,M,D,".grid.Z")
#   download.file(t ,temp, mode="wb")
#   file.copy(temp, "tmp.Z",overwrite=TRUE)
#   system("uncompress tmp.Z")
#   tmax <- raster("tmp")
#   
#   temp2 <- tempfile()#prep vapour pressure data
#   vp <- paste0("http://www.bom.gov.au/web03/ncc/www/awap/vprp/vprph15/daily/grid/0.05/history/nat/",Y,M,D,Y,M,D,".grid.Z")
#   download.file(vp ,temp2, mode="wb")
#   file.copy(temp2, "tmp2.Z",overwrite=TRUE)
#   system("uncompress tmp2.Z")
#   VP <- raster("tmp2")
#   
#   VPD<-VPDmaker(VP, tmax)
#   #writeRaster(VPD,paste0("~/Google Drive/PhD/Data/Spatial/Climate/VPD/VPD",Y,M,D,".asc"),overwrite=TRUE)
#   saveRDS(VPD, file = paste0("~/Google Drive/PhD/Data/Spatial/Climate/VPD2/VPD",Y,M,D,".rds"),
#     ascii = FALSE, compress = TRUE)
#   
#   
#   unlink(c("tmp","tmp2"))
#   rm("VP","tmax")
#   message(paste0(i," out of ",length(dates)))
# }
# #put the data together and find average VPD
# dat.dir<-"/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/VPD"
# files<-list.files(dat.dir,full.names = TRUE)
# 
# 
# vpd<-readRDS(files[1])
# for(i in 2:length(files)){
#   vpd<-vpd+readRDS(files[i])
#   message(i)
# }
# 
# vpd2<-vpd/length(files)

VPD<-raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/averageAnnualVPD-05degrees.asc", 
  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

VPDAlbers05<-projectRaster(VPD,km1,method='ngb')
VPD100<-aggregate(VPDAlbers05,100,'mean')

writeRaster(VPD100,paste0("/Users/daisy/Google Drive/PhD/ThermalStress/data/spatial/VPD100km",as.Date(Sys.time()),".asc"))
