remove_samples <- function(Dat,samples,droplevels = TRUE){
  if(class(Dat) != "Dataset")
    stop("ERROR: An object of class Dataset must be provided.",call.=TRUE)
  if(class(samples) != "character")
    stop("ERRPR: samples must be a character vector",call.=TRUE)
  if(!all(samples %in% colnames(Dat$Tab))){
    stop("ERROR: some of the elements in samples are not in the Dataset",call.=TRUE)
  }
  
  samples_to_keep <- !(colnames(Dat$Tab) %in% samples)
  
  tab <- Dat$Tab[ , samples_to_keep, drop = TRUE ]
  map <- Dat$Map[ samples_to_keep, , drop = TRUE ]
  tax <- Dat$Tax

  if(droplevels && !is.null(map)){
    map <- droplevels(map)
  }
  
  res <- create_dataset(Tab=tab,Map=map,Tax=tax) 
  return(res)
}
