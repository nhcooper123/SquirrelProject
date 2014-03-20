#Useful generic functions 

require(plyr)

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

#Making data subsets
get.subset <- function(data, subset.col, subset.no) { 
    data[which(data[[subset.col]]==unique(data[[subset.col]])[subset.no]),]
  }

#Functions to calculate PSR

#Calculating PSR for one host, with counting sp. where genus only occurs once
count.parasites <- function(x) {
  genus <- sub(" .*$", "", x)
  species <- sub("^.* ", "", x)
  ans <- tapply(species, genus, function(spp)
         if (length(spp) == 1) 1 else length(setdiff(spp, "sp.")))
  sum(ans)
}

#PSR for all host species in a dataset
PSR.all <- function(data, parasite.col, host.col) {
  aggregate(data[[parasite.col]] ~ data[[host.col]], data, count.parasites)
}

#Ouput for PSR calculations
PSR.output <- function(output, subset.name=NULL) {
  names(output) <- c("host", paste("PSR", subset.name, sep = ""))
  return(output)
}

#Overall PSR function
PSR <- function(data, parasite.col, host.col, subset.name=NULL) {
  psr <- PSR.all(data, parasite.col, host.col)
  PSR.output(psr, subset.name)
  } 
}