remove_taxons <- function(Dat,taxons){
  if(class(Dat) != "Dataset")
    stop("ERROR: An object of class Dataset must be provided.",call.=TRUE)
  if(class(taxons) != "character")
    stop("ERRPR: samples must be a character vector",call.=TRUE)
  if(!all(taxons %in% row.names(Dat$Tab))){
    stop("ERROR: some of the elements in samples are not in the Dataset",call.=TRUE)
  }
  
  taxons_to_keep <- !(row.names(Dat$Tab) %in% taxons)
  
  tab <- Dat$Tab[ taxons_to_keep, , drop = FALSE]
  map <- Dat$Map
  tax <- Dat$Tax[ taxons_to_keep, , drop = FALSE]
  
  res <- create_dataset(Tab=tab,Map=map,Tax=tax)
  return(res)
}