column.ID <- function(data, column.name) { 
  which(names(data)==column.name)
  }	

PSR.output <- function(output, subset.name=NULL) {
  output <- data.frame(output)
  output$host <- rownames(output)
  names(output)[1] <- paste("PSR", subset.name, sep = "")
  return(output)
}

PSR.all <- function(data, parasite.col, host.col) {
  host<-column.ID(data, host.col)
  parasite<-column.ID(data, parasite.col)
  psr <- tapply(data[,parasite], data[,host], FUN=length)
  PSR.output(psr)
} 

PSR.subset <- function(data, parasite.col, host.col, subset.col) {
  host <- column.ID(data,host.col)
  parasite <- column.ID(data, parasite.col)
  subset <- column.ID(data, subset.col)
  
  for(i in seq_along(unique(subset.col))) {
    psr <- tapply(data[,parasite][data[,subset]==unique(subset.col)[i]], data[,host][data[,subset]==unique(subset.col)[i]], FUN=length)
    PSR.output(psr, subset.name = unique(subset.col)[i])
  }
}

PSR <- function(data, parasite.col, host.col, subset.col=NULL){
  if(subset.col!= NULL) {
    psr <- PSR.subset(data, parasite.col, host.col, subset.col)
  } else {	
	psr <- PSR.all(data, parasite.col, host.col)
  }	
}






