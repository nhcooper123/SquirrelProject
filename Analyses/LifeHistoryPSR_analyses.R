library(ape)
source("PSR_functions.R")
source("Useful_functions.R")
source("SamplingEffort_function.R")
set.seed(1)

#-------------------------------------------------------------------
#Prepping data for PGLS analyses

#Currently this only keeps lines of the data
#which have full species names for squirrels, 
#where the squirrel is in the tree,
#and parasites with transmission modes and types
#-------------------------------------------------------------------

ds <- read.delim("squirrel.data.txt", header = TRUE)
mammal.tree <- read.nexus("mammalST_MSW05_best_chrono.tre")
PT <- read.delim("PanTHERIA05.txt")
#LH data
#citation count data

#Identify variables
host <- column.ID(ds, "HostCorrectedName")
parasite <- column.ID(ds, "ParasiteCorrectedName")
type <- column.ID(ds, "ParasiteType")
#transmode <- column.ID(psr.data, "ParasiteType")

#Tidy up the data 
ds <- remove.sp(ds, host)
ds <- remove.blanks(ds, host)
ds <- remove.blanks(ds, parasite)
ds <- remove.blanks(ds, type)
#ds <- remove.blanks(ds, transmode)
ds <- remove.incomplete.data(ds, c(host,parasite,type))

#Match tree and data
ds <- replace.spaces(ds, host)
ds <- remove.missing.species.tree(mammal.tree, ds, host)
squirrel.tree <- remove.missing.species.data (mammal.tree, ds, host)

#Estimate PSR for all species, parasite types, and transmission modes
psr.data <- unique.pairs(ds, parasite, host)
PSR <- psr(psr.data, parasite, host)
PSR.type <- psr(psr.data, parasite, host, type)
#PSR.transmode <- psr(ds, parasite, host, transmode)

#Estimate sampling effort for each host species
#TO DO#

#Merge the datasets together
PSR.complete <- merge(PSR, PSR.type, by.x = "host", by.y = "host")


squirrel.data <- merge(PSR, PSR.type, by.x = "host", by.y = "host") 

#-------------------------------------------------------------------
#PGLS analyses
#-------------------------------------------------------------------

squirrel <- comparative.data()

model1 <- pgls(PSR ~ ln(BodySize) + ln(citations), data = squirrel, lambda = 'ML')
summary(model1)
plot(model1)

plot(PSR ~ ln(BodySize), data = squirrel)
