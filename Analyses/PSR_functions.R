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

#Ouput for PSR calculations
psr.output <- function(output, subset.name=NULL) {
  names(output) <- c("host", paste("PSR", subset.name, sep = ""))
  return(output)
}

#PSR for all host species in a dataset
psr.all <- function(data, parasite.col, host.col) {
  ans <- aggregate(data[[parasite.col]] ~ data[[host.col]], data, count.parasites)
  psr.output(ans)
}

#PSR for subsets. Produces a big dataframe with results for each subset
#NA where there are no parasites of the subtype available

psr.subset <- function(data, parasite.col, host.col, subset.col) {
  result <- data.frame(host=NA) 
  
  for(i in seq_along(unique(data[[subset.col]]))) {
    data.subset <- get.subset(data, subset.col, i)
    ans <- psr.all(data.subset, parasite.col, host.col)
    ans <- psr.output(ans, subset.name=unique(data[[subset.col]])[i])

    result <- merge(ans, result, by = "host", all.x = TRUE, all.y = TRUE)
  }
  result <- result[!is.na(result$host),]
  return(result)
}

#Overall PSR function
psr <- function(data, parasite.col, host.col, subset.col=NULL) {
  if(is.null(subset.col)) {
    psr.all(data, parasite.col, host.col)
  } else {
    psr.subset(data, parasite.col, host.col, subset.col)
  }
}