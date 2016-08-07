clean <- function(Dat,verbose=TRUE){
  # Remove samples and OTUs with zero reads
  if(class(Dat) != "Dataset")
    stop("ERROR: An object of class Dataset must be provided.",call.=TRUE)
  
  tab <- Dat$Tab != 0
  
  samples_to_keep <-colSums(tab) > 0
  taxons_to_keep <- rowSums(tab) > 0
  
  tab <- Dat$Tab[ taxons_to_keep, samples_to_keep]
  map <- Dat$Map[ samples_to_keep, , drop = FALSE]
  tax <- Dat$Tax[ taxons_to_keep, , drop = FALSE]
  
  res <- create_dataset(Tab=tab,Map=map,Tax=tax)
  
  if(verbose){
    cat(sum(!samples_to_keep),"samples removed\n")
    cat(sum(!taxons_to_keep),"taxons removed\n")
  }
  
  return(res)
}
