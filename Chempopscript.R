#Import regional datasets containing species counts (Do for all 7 regions)
angliataxondata<-read.csv('Anglian_taxon_data_final_SPECIES.csv')
#merge regional taxon datasets
taxondataallregions<-merge(angliataxondata,midlandstaxondata,northeasttaxondata,northwesttaxondata,southerntaxondata,thamestaxondata,southwesttaxondata,by='REGION')
#Import regional datasets containing biotic index values (Do for all 7 regions)
angliasamplemaster<-read.csv('SamplemasterAnglia.csv') 
#merge regional samplemaster datasets into one country samplemaster data set
samplemasterallregions<-merge(angliasamplemaster,midlandssamplemaster,northeastsamplemaster,northwestsamplemaster,southernsamplmaster,thamessamplemaster,southwestsamplemaster,by='REGION')
#Import SITE ID guide
IDlist<-read.csv('IDlists.csv')
#Match Analysis Ids in taxon data sets to Site IDs in index value datasets
taxondataallregions$SITEID<-IDs$SITE_ID[match(samplemasterallregions$ANALYSIS_ID, IDs$ANALYSIS_ID)]
#Creata A Site-Species Matrix 
Matrixallregions<-with(taxondataallregions,tapply(SumOfTOTAL_ABUNDANCE,list(SITE_ID,FE_SPGRP_NAME),sum))
#Remove NAs from matrix
Matrixallregions[is.na(Matrixallregions)]<-0
#Import Phylogenetic Tree
install.packages('ape')
library(ape)
allregionsphylogeny<-read.tree(file='allregionsphylogeny.nwk')
#Create Comparative.comm object
install.packages('pez')
library(pez)
allregionscomparativeobject<-comparative.comm(allregionsphylogeny,allregionsMATRIX)
#Calculate Faith's PD values
pezshapeoutput<-pez.shape(allregionscomparativeobject)
#Calculate SESmpd and SESmntd values
pezdispersionoutput<-pez.dispersion(allregionscomparativeobject)
#select relevant columns from dispersion output
library(dplyr)
pezdispersionoutput<-pezdispersionoutput%>%select(SITE_ID,ses.mntd.mntd.mntd.obs.z,ses.mpd.mpd.obs.z)
#merge pez dispersion and pez shape data
outputallpez<-merge(pezshapeoutput,pezddispersionoutput, by='SITE_ID')
#calculate mean values for biotic indices for each site
outputpezmeans<-outputallpez%>%group_by(SITE_ID)%>%summarize(meanWHPT_ASPT=mean(WHPT_ASPT),meanDEHLI=mean(DEHLI),meanPSIAB=mean(PSI_AB))
#Linear models for biotic indices against pez metrics (example uses Faith's D and mean WHPT ASPT)
pdmodel1<-lm(formula=meanWHPT_ASPT~pd.pd, data=outputpezmeans)
summary(pd1model)
#plot linear regression 
library(ggplot2)
pdgraph1<-dggplot(outputpezmeans,aes(x=pd.pd, y=meanWHPT_ASPT))+geom_point()+labs(x='Faith D', y='Mean WHPT ASPT')+theme_bw()+  theme(axis.line = element_line(colour = "black"))
#Perform for all combinations of metrics and biotic indices
#Plot graphs together
install.packages('patchwork')
library(patchwork)
pdgraph1/(pdgraph2+pdgraph3)
#Perform and visualize Spearman's Rank Test 
spearmanplot<-outputpezmeans%>%select(meanWHPT_ASPT,meanPSI_AB,meanDEHLI,ses.mntd.mntd.obs.z,pd.pd,ses.mpd.mpd.obs.z,CONDUCTIVITY,ALKALINITY,WIDTH,DEPTH,SAND,SILT_CLAY,PEBBLES_GRAVEL,BOULDERS_COBBLES,ALTITUDE,SLOPE, DIST_FROM_SOURCE,DISCHARGE)
ggcorr(spearmanplot, method=c("pairwise", "spearman"), label=TRUE, label_size = 2, layout.exp = 10, hjust=1)
#Z-standardise physical variable values (Do for all 12 variables)
s.altitude<-scale(outputpezmeans$ALTITUDE) 
#Perform Akaike Information Criteria Model Averaging (repeat for each pd metric)
install.packages('MuMin')
library(MuMin)
pdmodeldata<- outputpezmeans%>%select(s.altitude,s.slope,s.distance,s.discharge,s.width,s.depth,s.boulders,s.pebbles,s.sand,s.silt,s.alkalinity+s.conductivity)
pdmodeldata<-na.omit(pdmodeldata)
pdmodel<-lm(pd.pd~s.altitude+s.slope+s.distance+s.discharge+s.width+s.depth+s.boulders+s.pebbles+s.sand+s.silt+s.alkalinity+s.conductivity, data=outputpezmeans, na.action='na.pass')
pdmodeldredged<-dredge(pdmodel)
pdmodelsummary<-summary(model.avg(dredgemodel2,subset=delta<4,fit=TRUE))





