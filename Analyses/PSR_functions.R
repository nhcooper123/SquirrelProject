#Useful generic functions 

#Identify column numbers
column.ID <- function(data, column.name) { 
  which(names(data)==column.name)
}	

#Functions for tidying up the data prior to PSR calculations:

#Remove replicated records of host-parasite pairs
unique.pairs <- function(data, parasite.col, host.col) {
  data[which(!duplicated(data[,c(host.col,parasite.col)])),]
}

#Remove all species with sp. not full Latin Binomial
remove.sp <- function(data, species.col) {
  data[-grep("sp.$", as.character(data[,species.col]),]
}

#Remove blank species names
remove.blanks <- function(data, species.col) {
  data <- data[-grep("sp.$", as.character(data[,species.col]),]###FIX####
}

#Functions to calculate PSR for data & subsets

#PSR for all species
PSR.all <- function(data, parasite.col, host.col) {
  genus <- sub(" .*$", "", x)
  species <- sub("^.* ", "", x)
  tapply(data[,parasite.col], data[,host.col], function(spp)
    if (length(spp) == 1) 1 else length(setdiff(spp, "sp.")))
}

  count.parasites <- function(x) {
  genus <- sub(" .*$", "", x)
  species <- sub("^.* ", "", x)
  tapply(species, genus, function(spp)
         if (length(spp) == 1) 1 else length(setdiff(spp, "sp.")))
} 

#PSR for just a subset of species
PSR.subset <- function(data, parasite.col, host.col, subset.col) {
  host <- column.ID(data,host.col)
  parasite <- column.ID(data, parasite.col)
  subset <- column.ID(data, subset.col)
  
  for(i in seq_along(unique(subset.col))) {
    psr <- tapply(data[,parasite][data[,subset]==unique(subset.col)[i]], data[,host][data[,subset]==unique(subset.col)[i]], FUN=length)
    PSR.output(psr, subset.name = unique(subset.col)[i])
  }
}



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








#To write
#Code to make output useable
#code that can deal with sp.s in Genus and Species names
#code that makes the dataset unique for pairs of species and genera
#Reminder to people to check spelling and taxonomy




