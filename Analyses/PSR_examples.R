#Examples for using PSR functions
#Data is available in SquirrelProject/Data

source("PSR_functions.R")

set.seed(1)
ds <- read.delim("squirrel.data.txt", header = TRUE)

#1: Clean up your dataset
#Identify your variables
#Remove replicated records of host-parasite pairs
#Remove host species with sp. not full Latin Binomial
#Remove blank species names

host <- column.ID(ds, "HostCorrectedName")
parasite <- column.ID(ds, "ParasiteCorrectedName")

ds <- unique.pairs(ds, parasite, host)
ds <- remove.sp(ds, host)
ds <- remove.blanks(ds, host)
ds <- remove.blanks(ds, parasite)

#2: PSR for all types of parasite

PSR <- psr(ds, parasite, host)
PSR
hist(PSR$PSR)

#3: PSR for different parasite types
#Identify column to subset data with
#Remove blank values for subset column

type <- column.ID(ds, "ParasiteType")
ds <- remove.blanks(ds, type)

PSR.type <- psr(ds, parasite, host, type)
PSR.type
names(PSR.type)
hist(PSR.type$PSRVirus)
