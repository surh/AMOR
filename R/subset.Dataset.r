#' Subset samples from a Dataset object.
#' 
#' Takes a Dataset object and a logical expression and returns a Dataset with only the samples
#' that meet the logical expression criteria
#' 
#' If the logical expression evaluation returns an NA value, it is considered equivalent to
#' FALSE
#' 
#' @param x A \code{Dataset} object.
#' @param subset Expression indicating samples to subset. See subset in \link{subset}
#' @param drop Logical indicating whether to drop variable levels that disappear.
#' @param clean Logical indicating whether to use the \link{clean} function on the resulting
#' dataset prior to returning it.
#' 
#' @author Sur from Dangl Lab
#' 
#' @return A Dataset object
subset.Dataset <- function(x,subset,drop = FALSE,clean = FALSE){
  
  e <- substitute(subset)
  r <- eval(e, x$Map, parent.frame())
  if (!is.logical(r)) 
    stop("'subset' must be logical")
  r <- r & !is.na(r) # removes NAs
  
  Map <- x$Map[ r, , drop = FALSE ]
  if(drop){
    Map <- droplevels(Map)
  }
  
  Tab <- x$Tab[ , row.names(Map), drop = FALSE ]
  
  Dat <- create_dataset(Tab = Tab, Map = Map, Tax = x$Tax)
  
  if(clean){
    Dat <- clean(Dat)
  }
  
  return(Dat)
}