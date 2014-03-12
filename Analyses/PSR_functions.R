#Identify column numbers
column.ID <- function(data, column.name) { 
  which(names(data)==column.name)
}	

#Remove replicated records of host-parasite pairs
unique.pairs <- function(data, parasite.col, host.col) {
  host <- column.ID(data,host.col)
  parasite <- column.ID(data, parasite.col)
  data <- data[which(duplicated(data[,c(host,parasite)])),]
  return(data)
}

#Remove all species with sp. not full Latin Binomial
remove.sp <- function(data, column.name) {
  column<-columnID(data,column.name)
  data[-grep("sp.$", as.character(data[,column.no]),]
}

#Remove blank species names
remove.blanks <- function(data, column.name) {
  data<-subset(data, column.name != "")
}

#Ouput for PSR calculations
PSR.output <- function(output, subset.name=NULL) {
  output <- data.frame(output)
  output$host <- rownames(output)
  names(output)[1] <- paste("PSR", subset.name, sep = "")
  return(output)
}

#PSR for all species
PSR.all <- function(data, parasite.col, host.col) {
  host<-column.ID(data, host.col)
  parasite<-column.ID(data, parasite.col)
  psr <- tapply(data[,parasite], data[,host], FUN=length)
  PSR.output(psr)
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




