require(caper)
source("PSR_functions.R")
source("Useful_functions.R")
source("SamplingEffort_functions.R")
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
LHdata <- read.delim("SquirrelLifeHistory.txt")
sampling.effort <- read.delim("squirrel.sampling.txt", header = TRUE)

#Identify variables
host <- column.ID(ds, "HostCorrectedName")
parasite <- column.ID(ds, "ParasiteCorrectedName")
type <- column.ID(ds, "ParasiteType")
vector <- column.ID(ds, "Vector")
intermediate <- column.ID(ds, "Intermediate")
nonclose <- column.ID(ds, "Nonclose")
vertical <- column.ID(ds, "Vertical")
sexual <- column.ID(ds, "Sexual")
close <- column.ID(ds, "CloseT")

#Tidy up the data 
ds <- remove.sp(ds, host)
ds <- remove.blanks(ds, host)
ds <- remove.blanks(ds, parasite)
ds <- remove.blanks(ds, type)
ds <- remove.incomplete.data(ds, c(host,parasite,type, 
	                               vector, intermediate, nonclose,
	                               vertical, sexual, close))

#Match tree and data
ds <- replace.spaces(ds, host)
ds <- remove.missing.species.tree(mammal.tree, ds, host)
squirrel.tree <- remove.missing.species.data(mammal.tree, ds, host)

#Estimate PSR for all species, parasite types, and transmission modes
#Then merge to make a dataset with all PSR results
psr.data <- unique.pairs(ds, parasite, host)
PSR <- psr(psr.data, parasite, host)
PSR.type <- psr(psr.data, parasite, host, type)
PSR.vector <- psr(ds, parasite, host, vector, binary=TRUE)
PSR.vertical <- psr(ds, parasite, host, vertical, binary=TRUE)
PSR.intermediate <- psr(ds, parasite, host, intermediate, binary=TRUE)
PSR.close <- psr(ds, parasite, host, close, binary=TRUE)
PSR.nonclose <- psr(ds, parasite, host, nonclose, binary=TRUE)
PSR.sexual <- psr(ds, parasite, host, sexual, binary=TRUE)

list.of.PSR.results <- list(PSR, PSR.type, PSR.vector, PSR.vertical,
							PSR.intermediate, PSR.sexual, PSR.close,
							PSR.nonclose)

PSR.complete <- Reduce(function(...) merge(..., by="host", all=TRUE), 
	                   list.of.PSR.results)
PSR.complete[is.na(PSR.complete)]<-0

#Estimate sampling effort for each host species
hostname <- column.ID(sampling.effort, "HostCorrectedName")
samples <- sampling.occ(sampling.effort, hostname)

#Merge the datasets together
squirrel.data <- merge(PSR.complete, samples, by="host", all=TRUE) 

squirrel.data <- merge(squirrel.data, LHdata, by.x="host", 
	                   by.y="MSW05_Binomial" all.x=TRUE) 

#And finally...make comparative data object for PGLS

squirrel <- comparative.data(phy=squirrel.tree, data=squirrel.data, 
	                         names.col=squirrel.data[[host]])

#-------------------------------------------------------------------
#PGLS analyses
#-------------------------------------------------------------------
#Example: Total parasite species richness and body size
#With sampling occassions to control for sampling effort

model1 <- pgls(log(PSR) ~ log(BodyMass_g) + log(sampling.occ), data=squirrel, lambda='ML')
summary(model1)
plot(model1)

plot(log(PSR) ~ log(BodyMass_g), data=squirrel)

#Example: Helminth parasite species richness and body size
#With sampling occassions to control for sampling effort

model2 <- pgls(log(PSRhelminth) ~ log(BodyMass_g) + log(sampling.occ), data=squirrel, lambda='ML')
summary(model2)
plot(model2)
