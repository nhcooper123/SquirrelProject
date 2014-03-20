#Useful generic functions 

#Identify column numbers
column.ID <- function(data, column.name) { 
  which(names(data)==column.name)
}	

#Functions for tidying up the data prior to PSR calculations:

#Remove replicated records of host-parasite pairs
unique.pairs <- function(data, parasite.col, host.col) {
  data[which(!duplicated(data[c(host.col,parasite.col)])),]
}

#Remove all species with sp. not full Latin Binomial
remove.sp <- function(data, species.col) {
  data[-grep("sp.$", as.character(data[[species.col]])),]
}

#Remove blank species names
remove.blanks <- function(data, species.col) {
  data[which(data[species.col] != ""),]
}

#Functions to calculate PSR for data & subsets

#Calculating PSR for one host, with counting sp. where genus only occurs once
count.parasites <- function(parasite) {
  genus <- sub(" .*$", "", parasite)
  species <- sub("^.* ", "", parasite)
  sum(tapply(species, genus, function(spp)
         if (length(spp) == 1) 1 else length(setdiff(spp, "sp."))))
}

#PSR for all host species
PSR.all <- function(data, parasite.col, host.col) {
  tapply(data[[parasite.col]], data[[host.col]], count.parasites)
}

#PSR for just a subset of species
PSR.subset <- function(data, parasite.col, host.col, subset.col) { 
  for(i in seq_along(unique(data[[subset.col]]))) {
    data.subset <- data[which(data[[subset.col]]==unique(data[[subset.col]])[i]),]
    print(PSR.all(data.subset, parasite.col, host.col))
  }
}

SPELT.data[i,3] is slower than SPELT.data[[3]][i]

#Ouput for PSR calculations
PSR.output <- function(output, subset.name=NULL) {
  output <- data.frame(output)
  output$host <- rownames(output)
  names(output)[1] <- paste("PSR", subset.name, sep = "")
  return(output)
}

#Combined function to run any PSR calculation
PSR <- function(data, parasite.col, host.col, subset.col=NULL) {
  if(subset.col!= NULL) {
    psr <- PSR.subset(data, parasite.col, host.col, subset.col)
  } else {	
	psr <- PSR.all(data, parasite.col, host.col)
  }	
}







psr <- tapply(data[,parasite], data[,host], FUN=length)
#To write
#Code to make output useable
#code that can deal with sp.s in Genus and Species names
#code that makes the dataset unique for pairs of species and genera
#Reminder to people to check spelling and taxonomy




