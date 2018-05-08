#' Normalize
#' 
#' Takes a dataset object and normalizes the abundance
#' matrix by a variable from Map or a vector provided.
#' 
#' @param Dat A Dataset object.
#' @param norm Either a string indicating the variable name from
#' Map to use for normalization, or a numeric vector with normalization
#' values
#' 
#' @export
normalize <- function(Dat, norm = NULL){
  if(class(Dat) != "Dataset")
    stop("ERROR: Dat must be of class Dataset")
  
  Tab <- Dat$Tab
  if(length(norm) == 1 && is.character(norm)){
    norm <- Dat$Map[,norm]
  }else if(length(norm) == ncol(Tab) && is.numeric(norm) && all(norm > 0)){
    # Proceed
  } else {
    stop("ERROR: norm must be a positive numeric vector of length equal to sample number,\nor a character string indicating the variable name to use as normalization factor")
  }
  
  Tab <- t(100*t(Tab) / norm)
  Dat <- create_dataset(Tab = Tab, Map = Dat$Map, Tax = Dat$Tax)
  return(Dat)
}