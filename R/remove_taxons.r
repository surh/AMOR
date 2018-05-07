#' Remove taxa from a Dataset object.
#' 
#' Remove taxa from a Dataset object.
#'
#' @param Dat A Dataset object.
#' @param taxons String vector indicating
#' which taxa to remove. Taxon IDs must be
#' equal to the row names of the abundance
#' table on Dat.
#'
#' @return A Dataset object.
#' @export
#' @author Sur Herrera Paredes
#' @seealso \link{remove_samples}
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' Dat <- create_dataset(Rhizo, Rhizo.map)
#' Dat <- remove_samples(Dat,
#'                       colnames(Dat$Tab)[1:50])
#' Dat <- remove_taxons(Dat,
#'                      row.names(Dat$Tab)[1:115])
#' sum(colSums(Dat$Tab) == 0)
#' sum(rowSums(Dat$Tab) == 0)
#' 
#' Dat$Tab
#' Dat.clean <- clean(Dat)
#' Dat.clean$Tab
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