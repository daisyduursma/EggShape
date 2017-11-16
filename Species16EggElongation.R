

# test to see if elongation changes within a species due to differences in VPD
#
# **Removed all values from QVM - the were incorrectly entered into the .csv file
#
# mixed effect model with VPD as fixed effect and nest ID as random effects
#
# "Linear mixed-effects models were used to assess how within species egg elongation changed in response to VPD. ANOVA was calculated using an Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger approximation for degrees-of-freedom). In our models for each species  VPD was a fixed effect and nest ID was a random effect to account for differences in egg shape within a single nest."
#

rm(list = ls())
library(lme4)
library(raster)
library(car)


#Data
spDat <- read.csv("~/Google Drive/PhD/ThermalStress/data/151124_otherSpecies_egg_collection_details_rearranged.csv")
VPDann<-raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/averageAnnualVPD.asc")

#species of interest
species<-c("Cracticus tibicen","Petrochelidon ariel","Zosterops lateralis",
           "Anthus australis","Rhipidura leucophrys","Grallina cyanoleuca",
           "Petrochelidon nigricans","Taeniopygia bichenovii","Pachycephala rufiventris",
           "Rhipidura albiscapa","Pomatostomus temporalis","Melanodryas cucullata",
           "Manorina melanocephala","Anthus novaeseelandiae")

#rename
nms<-c(names(spDat)[1:12],paste(rep(c("egg_area","egg_length","egg_width"),11),rep(1:11,each=3),sep="."),names(spDat)[46])
colnames(spDat)<-nms

#change to wide format
longDat<-list()
for (l in 1:11){
  df1<-spDat[,c(names(spDat)[1:12],paste(c("egg_area","egg_length","egg_width"),rep(l,3),sep="."),names(spDat)[46])]
  colnames(df1)<-c(names(spDat)[1:12],"egg_area","egg_length","egg_width",names(spDat)[46])
  longDat[[l]]<-df1
}
longDat<-as.data.frame(do.call(rbind,longDat))
longDat<-subset(longDat,!is.na(longDat$egg_area))#remove NA values
longDat<-subset(longDat,!is.na(longDat$egg_length))#remove NA values

#Extract VPD at the locations
longDat$VPD<-extract(VPDann,cbind(longDat$decimalLongitude,longDat$decimalLatitude))

#for cells falling just off the coast use bilinear interpolation to get VPD values
longDatNA<-subset(longDat,is.na(longDat$VPD))
longDat<-subset(longDat,!is.na(longDat$VPD))
longDatNA$VPD<-extract(VPDann,cbind(longDatNA$decimalLongitude,longDatNA$decimalLatitude),method='bilinear')
longDatNA<-subset(longDatNA,!is.na(longDatNA$VPD))
longDat<-rbind(longDat,longDatNA)
longDat$elongation<-with(longDat,egg_length/egg_width)


#loop through each species and find out if egg elongation changes with VPD,
#we expect some species will have decreased elongation with increasing VPD
#some nests have multiple eggs so treat nest ID as random effect

modelDat<-list()
for (sp in 1:length(species)){
  df2<-subset(longDat, Species == species[sp])
  observationCount<-nrow(df2)
  uniqueNests<-length(unique(df2$id))
  #mixed effect model to control for multiple eggs in a nest
  mod1<-lmer(elongation ~ VPD + (1+VPD|id),data=df2,REML=TRUE)#add nest ID
    # #check on spatial autocorrelation
      # C_lmEl <- correlog(df2$decimalLongitude, df2$decimalLatitude, residuals(modBP1),
      #                    na.rm=T, increment=100, resamp=0,latlon=TRUE)
      #
      # # plot(0,type="n",col="black",ylab="Moran's I",xlab="lag distance",xlim=c(0,4600),ylim=c(-1,1),main = "Elongation")
      # abline(h=0,lty="dotted")
      # lines(C_lmEl$correlation~C_lmEl$mean.of.class,col="red",lwd=2)
      # plot(VPDann)
      # points(df2$decimalLongitude, df2$decimalLatitude, col=c("blue",
      #                          "red")[sign(resid(modBP1))/2+1.5], pch=19,
      #      cex=abs(resid(modBP1))/max(resid(modBP1))*2, add=TRUE)

  # qqPlot(residuals(modBP1))
  #text(120,-10,paste(sp))
  p<-round(Anova(mod1)$"Pr(>Chisq)",digits = 3)#get p values from lmer
    #text(130,-40,p)
  modAnova<-Anova(mod1, test.statistic="F")
  f<-modAnova$F
  df<-modAnova$Df
  dfRes<-modAnova$"Df.res"

  modelDat[[sp]]<-cbind(species[sp],observationCount,uniqueNests,p,f,df,dfRes)
}

do.call("rbind",modelDat)

# If I wanted to report results for sp <- 14

#F(1, 4.23) = 0.5495518, p = 0.37
