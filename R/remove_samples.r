#' Remove samples from a Dataset object.
#' 
#' Remove samples from a Dataset object.
#'
#' @param Dat A Dataset object
#' @param samples String vector indicating which
#' samples to remove. Sample IDs must be equal to
#' the column headers of the abundance table on Dat.
#' @param droplevels Logical indicating whether
#' levels that desappear from the mapping file should
#' be removed.
#'
#' @return A Dataset object
#' @export
#' @author Sur Herrera Paredes
#' @seealso \link{remove_taxons}
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' Dat <- create_dataset(Rhizo,Rhizo.map)
#' Dat <- remove_samples(Dat = Dat,
#'                       samples = colnames(Dat$Tab)[1:50])
#' Dat <- remove_taxons(Dat = Dat,
#'                      taxons = row.names(Dat$Tab)[1:115])
#' 
#' sum(colSums(Dat$Tab) == 0)
#' sum(rowSums(Dat$Tab) == 0)
#' 
#' Dat$Tab
#' Dat.clean <- clean(Dat)
#' Dat.clean$Tab
remove_samples <- function(Dat, samples, droplevels = TRUE){
  if(class(Dat) != "Dataset")
    stop("ERROR: An object of class Dataset must be provided.",call.=TRUE)
  if(class(samples) != "character")
    stop("ERRPR: samples must be a character vector",call.=TRUE)
  if(!all(samples %in% colnames(Dat$Tab))){
    stop("ERROR: some of the elements in samples are not in the Dataset",call.=TRUE)
  }
  
  samples_to_keep <- !(colnames(Dat$Tab) %in% samples)
  
  tab <- Dat$Tab[ , samples_to_keep, drop = FALSE ]
  map <- Dat$Map[ samples_to_keep, , drop = FALSE ]
  tax <- Dat$Tax

  if(droplevels && !is.null(map)){
    map <- droplevels(map)
  }
  
  res <- create_dataset(Tab=tab,Map=map,Tax=tax) 
  return(res)
}