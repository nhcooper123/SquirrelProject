# Useful Generic Functions

# Identify column numbers
column.ID <- function(data, column.name) { 
  which(names(data) == column.name)
} 

# Remove blanks from a column that should be complete
remove.blanks <- function(data, complete.col) {
  data[which(data[complete.col] != ""), ]
}

# Making data subsets
get.subset <- function(data, subset.col, subset.no) { 
    data[which(data[[subset.col]] == unique(data[[subset.col]])[subset.no]),]
}

# Making data subsets for binary categories
get.subset.binary <- function(data, subset.col) { 
    data[which(data[[subset.col]] == 1), ]
}

# Remove incomplete data for specific variables
remove.incomplete.data <- function(data, variable.col.list) {
  id <- complete.cases(data[, c(variable.col.list)])
  data <- data[id, ]
  return(data)
}

# Replace spaces in speciesnames with _
replace.spaces <- function(data, speciesnames.col) {
  data[[speciesnames.col]] <- gsub(" ", "_", data[[speciesnames.col]])
  return(data)
}

# ID species in tree that are not in the data
id.missing.tree <- function(phy, data, speciesnames.col) {
  setdiff(phy$tip.label, data[, speciesnames.col])
}

# ID species in data that are not in the tree
id.missing.data <- function(phy, data, speciesnames.col) {
  setdiff(data[, speciesnames.col], phy$tip.label)
}    

# Remove missing species from tree
remove.missing.species.tree <- function(phy, data, speciesnames.col) {
  tree.not.data <- id.missing.tree(phy, data, speciesnames.col)
  if (length(tree.not.data) > 0) {
    phy <- drop.tip(phy, tree.not.data)
  }
  return(phy)
}

# Remove missing species from data
remove.missing.species.data <- function(phy, data, speciesnames.col) {
  data.not.tree <- id.missing.data(phy, data, speciesnames.col)
  if (length(data.not.tree) > 0) {
    matches <- match(data[, speciesnames.col], data.not.tree, nomatch = 0)
    data <- subset(data, matches == 0)
  }
  return(data)
}