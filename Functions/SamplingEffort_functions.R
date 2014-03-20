#Sampling effort for parasite analyses

sampling.occ <- function(data, host.col) {
  data$ID <- 1:nrow(data)
  ans <- aggregate(data$ID ~ data[[host.col]], data, length)
  names(ans) <- c("host", "sampling.occasions")
  return(ans)
}