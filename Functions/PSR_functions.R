#Functions for tidying up the data prior to PSR calculations:

#Remove replicated records of host-parasite pairs
unique.pairs <- function(data, parasite.col, host.col) {
  data[which(!duplicated(data[c(host.col,parasite.col)])),]
}

#Remove all species with sp. not full Latin Binomial
remove.sp <- function(data, species.col) {
  if(length(grep("sp.$", as.character(data[[species.col]])))>0) {
    data[-grep("sp.$", as.character(data[[species.col]])),]
  }
  return(data)
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
#Adds zeros where there are no parasites of the subtype available

psr.subset <- function(data, parasite.col, host.col, subset.col) {
  result <- data.frame(host=NA) 
  
  for(i in seq_along(unique(data[[subset.col]]))) {
    data.subset <- get.subset(data, subset.col, i)
    ans <- psr.all(data.subset, parasite.col, host.col)
    ans <- psr.output(ans, subset.name=unique(data[[subset.col]])[i])

    result <- merge(ans, result, by = "host", all.x = TRUE, all.y = TRUE)
  }
  result <- result[!is.na(result$host),]
  result[is.na(result)]<-0
  return(result)
}

psr.subset.binary <- function(data, parasite.col, host.col, subset.col) {
  data.subset <- get.subset.binary(data, subset.col)
  ans <- psr.all(data.subset, parasite.col, host.col)
  result <- psr.output(ans, subset.name=names(data[subset.col]))
  return(result)
}

#Overall PSR function
psr <- function(data, parasite.col, host.col, subset.col=NULL, binary=NULL) {
  if(is.null(subset.col)) {
    psr.all(data, parasite.col, host.col)
  } else {
    if(is.null(binary)) {
      psr.subset(data, parasite.col, host.col, subset.col)
    } else {
      psr.subset.binary(data, parasite.col, host.col, subset.col)
    }
  }
}