column.ID <- function(data, column.name) { 
  which(names(data)==substitute(column.name))
  }	

PSR.output <- function(output, subset.name=NULL) {
  output <- data.frame(output)
  output$host <- rownames(output)
  names(output)[1] <- paste("PSR", subset.name, sep = "")
  return(output)
}