# Running analyses of PSR ~ life history variables
# Currently this only keeps lines of the data
# which have full species names for squirrels, 
# where the squirrel is in the tree,
# and parasites with transmission modes and types

require(caper)
source("PSR_functions.R")
source("Useful_functions.R")
source("SamplingEffort_functions.R")

# Read in datasets (will need to setwd or similar)

ds <- read.delim("squirrel.data.txt", header = TRUE)
american.species <- read.delim("SquirrelsofAmerica.txt", header = FALSE)
mammal.tree <- read.nexus("mammalST_MSW05_best_chrono.tre")
LHdata <- read.delim("SquirrelLifeHistoryData.txt", header = TRUE)
sampling.effort <- read.delim("SquirrelSamplingEffort.txt", header = TRUE)

# 1: FULL DATASET cleaning
# Identify variables
host <- column.ID(ds, "HostCorrectedName")
parasite <- column.ID(ds, "ParasiteCorrectedName")
type <- column.ID(ds, "ParasiteType")
vector <- column.ID(ds, "Vector")
intermediate <- column.ID(ds, "Intermediate")
nonclose <- column.ID(ds, "Nonclose")
vertical <- column.ID(ds, "Vertical")
sexual <- column.ID(ds, "Sexual")
close <- column.ID(ds, "CloseT")

# Subset out just the american species
american.sp <- match(ds[[host]], american.species[[1]], nomatch = 0)
ds <- subset(ds, american.sp != 0)

# Tidy up the data 
ds <- remove.sp(ds, host)
ds <- remove.blanks(ds, host)
ds <- remove.blanks(ds, parasite)
ds <- remove.blanks(ds, type)
ds <- remove.incomplete.data(ds, c(host,parasite,type, 
	                               vector, intermediate, nonclose,
	                               vertical, sexual, close))

# Match tree and data
ds <- replace.spaces(ds, host)
ds <- remove.missing.species.data(mammal.tree, ds, host)
squirrel.tree <- remove.missing.species.tree(mammal.tree, ds, host)

# 2: Extract PSR DATA
# Estimate PSR for all species, parasite types, and transmission modes
# Then merge to make a dataset with all PSR results
psr.data <- unique.pairs(ds, parasite, host)
PSR <- psr(psr.data, parasite, host)
PSR.type <- psr(psr.data, parasite, host, type)
PSR.vector <- psr(psr.data, parasite, host, vector, binary = TRUE)
PSR.vertical <- psr(psr.data, parasite, host, vertical, binary = TRUE)
PSR.intermediate <- psr(psr.data, parasite, host, intermediate, binary = TRUE)
PSR.close <- psr(psr.data, parasite, host, close, binary = TRUE)
PSR.nonclose <- psr(psr.data, parasite, host, nonclose, binary = TRUE)
PSR.sexual <- psr(psr.data, parasite, host, sexual, binary = TRUE)

list.of.PSR.results <- list(PSR, PSR.type, PSR.vector, PSR.vertical,
							PSR.intermediate, PSR.sexual, PSR.close,
							PSR.nonclose)

PSR.complete <- Reduce(function(...) merge(..., by = "host", all = TRUE), 
	                   list.of.PSR.results)
PSR.complete[is.na(PSR.complete)] <- 0

# 3: SAMPLING EFFORT DATA
# Identify variables
hostname <- column.ID(sampling.effort, "HostCorrectedName")

# Estimate sampling effort for each host species
samples <- sampling.occ(sampling.effort, hostname)
samples <- replace.spaces(samples, hostname)

# 4: Life History data
# Identify variables and replace spaces in species names
speciesname <- column.ID(LHdata, "MSW05_Binomial")
LHdata <- replace.spaces(LHdata, speciesname)

# 5: COMBINE DATASETS
# Merge the datasets together
squirrel.data <- merge(PSR.complete, samples, by = "host", all.x = TRUE) 

squirrel.data <- merge(squirrel.data, LHdata, by.x = "host", 
	                   by.y = "MSW05_Binomial", all.x = TRUE) 

# 5: And finally...make comparative data object for PGLS

squirrel <- comparative.data(phy = squirrel.tree, data = squirrel.data, 
	                        names.col = host, na.omit = FALSE)

#-------------------------------------------------------------------
# PGLS analyses
#-------------------------------------------------------------------
# Example: Total parasite species richness and body size
# With sampling occasions to control for sampling effort

model1 <- pgls(log(PSR) ~ log(AdultBodyMass_g) + log(sampling.occasions), 
	           data = squirrel, lambda = 'ML')
summary(model1)
plot(model1)

plot(log(PSR) ~ log(AdultBodyMass_g), data = squirrel$data)

# Example: Helminth parasite species richness and body size
# With sampling occasions to control for sampling effort
# Note that PSRHelminth can = 0 so we add 0.1 to log it

model2 <- pgls(log(PSRHelminth + 0.1) ~ log(AdultBodyMass_g) + log(sampling.occasions), 
	           data = squirrel, lambda = 'ML')
summary(model2)

# Example: Total parasite species richness and coloniality
# With sampling occasions to control for sampling effort

model3 <- pgls(log(PSR) ~ Colonial + log(sampling.occasions), 
	           data = squirrel, lambda = 'ML')
summary(model3)
