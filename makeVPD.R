library(raster)
dat.dir<-"/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/VPD"
files<-list.files(dat.dir,full.names = TRUE)


vpd<-readRDS(files[1])
for(i in 2:length(files)){
  vpd<-vpd+readRDS(files[i])
  message(i)
}

vpd2<-vpd/length(files)
writeRaster(vpd2,"/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/averageAnnualVPD-05degrees.asc")
