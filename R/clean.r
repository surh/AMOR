#' Remove samples and taxons from a Dataset
#' 
#' Functions that will take a Dataset object, and
#' remove samples and/or taxons with different criteria
#' 
#' clean will remove any sample (column) or taxon (row) fromt your Dataset,
#' that has no observations. Some functions like PCA and some linear models 
#' won't work if one of the rows or columns of the matrix are all zero,
#' and so it is always a good idea to "clean" your dataset.
#' 
#' remove_samples and remove_taxons take respectively a list of samples or
#' taxons to remove from the Dataset.
#'
#' @param Dat A Dataset object
#' @param verbose logical indicating whether to print info about the results.
#'
#' @return A Dataset object
#' @author Sur Herrera Paredes
#' @seealso \link{remove_taxons}, \link{remove_samples}
#' @export
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' Dat <- create_dataset(Rhizo,Rhizo.map)
#' Dat <- remove_samples(Dat,colnames(Dat$Tab)[1:50])
#' Dat <- remove_taxons(Dat,row.names(Dat$Tab)[1:115])
#' sum(colSums(Dat$Tab) == 0)
#' sum(rowSums(Dat$Tab) == 0)
#' Dat$Tab
#' Dat.clean <- clean(Dat)
#' Dat.clean$Tab
clean <- function(Dat,verbose=TRUE){
  # Remove samples and OTUs with zero reads
  if(class(Dat) != "Dataset")
    stop("ERROR: An object of class Dataset must be provided.",call.=TRUE)
  
  tab <- Dat$Tab != 0
  
  samples_to_keep <-colSums(tab) > 0
  taxons_to_keep <- rowSums(tab) > 0
  
  tab <- Dat$Tab[ taxons_to_keep, samples_to_keep, drop = FALSE]
  map <- Dat$Map[ samples_to_keep, , drop = FALSE]
  tax <- Dat$Tax[ taxons_to_keep, , drop = FALSE]
  
  res <- create_dataset(Tab=tab,Map=map,Tax=tax)
  
  if(verbose){
    cat(sum(!samples_to_keep),"samples removed\n")
    cat(sum(!taxons_to_keep),"taxons removed\n")
  }
  
  return(res)
}
