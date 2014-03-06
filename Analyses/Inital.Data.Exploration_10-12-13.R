#-------------------------------------
#Initial exploration of squirrel data
#10-12-13
#Natalie Cooper
#-------------------------------------

#------------------------------------------------------------
#Load required packages
#If you get errors you may need to install the packages first
#For each one:
#install.packages("maptools")
#------------------------------------------------------------

library(maptools)
gpclibPermit()
library(rgdal)
library(sp)
library(raster)
library(nlme)

library(PBSmapping)
library(rworldmap)
library(RColorBrewer)
library(rgeos)

library(ape)

#-----------------------------------
#Set directory you put the data into
#-----------------------------------

setwd("/Users/nataliecooper/Desktop/SquirrelProject")

ds<-read.delim("squirrel.data.txt")

#hello thomas!
#Hello again

str(ds)
names(ds)
head(ds)

#------------------------------------
#Some simple data cleaning
#------------------------------------

#Remove all host species with sp.
ds<-ds[-grep("sp.$", as.character(ds$HostCorrectedName)),]

#Remove blank host and parasite names
ds<-subset(ds, ParasiteCorrectedName != "")
ds<-subset(ds, HostCorrectedName != "")

#-------------------------------------
#Summary data
#-------------------------------------

length(unique(ds$HostCorrectedName))

length(unique(ds$ParasiteCorrectedName))

#--------------------------------------------
#Parasite species richness (psr) calculations
#Note that we're only interested in unique
#combinations of parasite + host
#So we need to cut out extra data and select
#unique entries only first
#-------------------------------------------

#Remove locality data (columns 3 through 7)
ds.psr<-ds[,-c(3:7)]

#Select unique combinations
ds.psr<-unique(ds.psr)

#Now calculate PSR
PSR_total<-with(ds.psr, aggregate(ParasiteCorrectedName, by = list(HostCorrectedName), FUN = length))

#To view results:
PSR_total

#Or in order of PSR:
PSR_total[rev(order(PSR_total$x)),]

#histogram of PSR values
hist(PSR_total$x)

#Which species has the most parasites?
PSR_total$Group.1[which(PSR_total$x == max(PSR_total$x))]

#We can also estimate PSR separately for different types of parasite...
PSR_arthropod<-with(ds.psr, aggregate(ParasiteCorrectedName[ParasiteType == "Arthropod"], by = list(HostCorrectedName[ParasiteType == "Arthropod"]),FUN = length))
PSR_bacteria<-with(ds.psr, aggregate(ParasiteCorrectedName[ParasiteType == "Bacteria"], by = list(HostCorrectedName[ParasiteType == "Bacteria"]),FUN = length))
PSR_helminth<-with(ds.psr, aggregate(ParasiteCorrectedName[ParasiteType == "Helminth"], by = list(HostCorrectedName[ParasiteType == "Helminth"]),FUN = length))
PSR_protozoa<-with(ds.psr, aggregate(ParasiteCorrectedName[ParasiteType == "Protozoa"], by = list(HostCorrectedName[ParasiteType == "Protozoa"]),FUN = length))
PSR_virus<-with(ds.psr, aggregate(ParasiteCorrectedName[ParasiteType == "Virus"], by = list(HostCorrectedName[ParasiteType == "Virus"]),FUN = length))

#or transmission modes...
PSR_close<-with(ds.psr, aggregate(ParasiteCorrectedName[CloseT == 1], by = list(HostCorrectedName[CloseT == 1]),FUN = length))
PSR_sexual<-with(ds.psr, aggregate(ParasiteCorrectedName[Sexual == 1], by = list(HostCorrectedName[Sexual == 1]),FUN = length))
PSR_vector<-with(ds.psr, aggregate(ParasiteCorrectedName[Vector == 1], by = list(HostCorrectedName[Vector == 1]),FUN = length))
PSR_intermediate<-with(ds.psr, aggregate(ParasiteCorrectedName[Intermediate == 1], by = list(HostCorrectedName[Intermediate == 1]),FUN = length))
PSR_nonclose<-with(ds.psr, aggregate(ParasiteCorrectedName[Nonclose == 1], by = list(HostCorrectedName[Nonclose == 1]),FUN = length))
PSR_vertical<-with(ds.psr, aggregate(ParasiteCorrectedName[Vertical == 1], by = list(HostCorrectedName[Vertical == 1]),FUN = length))

#------------------------------------
#Geographic distributions
#------------------------------------
#First create a new clean dataset
#with no missing coordinate data
#and data with high accuracy
#------------------------------------

#Remove data with no coordinate data
idx<-complete.cases(ds[,4:5])
ds.geographic<-ds[idx,]

#Remove inaccurate coordinates (Accuracy = 3 or 4 or Extent = 4)
#Can play around with options here
ds.geographic<-subset(ds.geographic, Accuracy < 3)
ds.geographic<-subset(ds.geographic, Extent_km2 < 4)

#-----------------------------------------------------
#convert to R Data Structure 'SpatialPointsDataFrame'
#-----------------------------------------------------

#create spatial points dataframe (SPDF)
#first thing needs to be long then lat, then the second is all the data
dsSPDF<-SpatialPointsDataFrame(ds.geographic[,4:5],data.frame(ds.geographic[,1:14]))

#make sure the projection is lat long
proj4string(dsSPDF)<-CRS("+proj=longlat")

#-------------------------
#Plot points on world map
#-------------------------
#Add world map data
data(wrld_simpl)

plot(wrld_simpl)
points(dsSPDF, col = "red", cex = 0.5, pch = 16)

#-----------------------------------
#Phylogenetic distribution
#----------------------------------
#Read in mammal phylogeny
#We will cut out non-squirrels next 
#-----------------------------------

mammal.tree<-read.nexus("mammalST_MSW05_best_chrono.tre")

#---------------------------------------------------------
#Add PanTHERIA ecological trait data
#This is an easy way of getting a list of squirrel species
#---------------------------------------------------------

PT<-read.delim("PanTHERIA05.txt")

#Select only species in the family Sciuridae
PT.squirrel<-subset(PT, MSW05_Family == "Sciuridae")

#Put underscores in spaces to match species names to those in tree
PT.squirrel$MSW05_Binomial<-gsub(" ", "_", PT.squirrel$MSW05_Binomial)

#-----------------------------------------------
#Remove non-squirrels from phylogeny
#-----------------------------------------------

squirrel.tree<-drop.tip(mammal.tree, setdiff(mammal.tree$tip.label, PT.squirrel$MSW05_Binomial))

#-----------------------------------------------
#Plot phylogeny for all squirrels
#Then add colours to squirrels we have data for
#-----------------------------------------------

plot(squirrel.tree, type = "fan", cex = 0.7)

#Put underscores in spaces to match species names to those in tree
ds$HostCorrectedName<-gsub(" ", "_", ds$HostCorrectedName)

#Make a new dataframe with a list of colours
#blue for species in the tree and our dataset
#grey for those not in our dataset
color.dataframe<-data.frame(species = squirrel.tree$tip.label, 
								indataset = match(squirrel.tree$tip.label, ds$HostCorrectedName, nomatch = 0), 
									colour.to.plot = rep(0, length(squirrel.tree$tip.label)))

color.dataframe$colour.to.plot[which(color.dataframe$indataset > 0)]<-"blue1"
color.dataframe$colour.to.plot[which(color.dataframe$indataset == 0)]<-"darkgray"
#See colors() for full list of colours in R

#plot tree
plot(squirrel.tree, type = "fan", cex = 0.7, tip.color = color.dataframe$colour.to.plot)

