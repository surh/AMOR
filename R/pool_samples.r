#' Pool samples
#' 
#' Pools abundances from samples according to a grouping factor.
#' 
#' 
#' Wrapper for \code{collapse_matrix}. This function is useful to calculate
#' per-group summary statistics per taxon.
#' 
#' The default method takes an abundance matrix and a grouping factor,
#' then applies the aggregating function FUN to the groups of samples
#' defined by the grouping factor.
#' 
#' The Dataset method takes a Dataset object and obtains the grouping factor
#' from the Map element.
#' 
#' @param x Either a numerical matrix or a Dataset object. For a numerical
#' matrix, samples are given as columns and taxa as rows.
#' @param groups For the default method this is a vector or factor 
#' specifying to which group each sample belongs. Vectors will be converted to
#' a factors with the \code{factor} function. It must be of the same length
#' as \code{ncol(x)}. For the Dataset method, this is a single character string or
#' numeric value indicating the column from x$Map to be used as a grouping factor.
#' @param FUN Function to apply when collapsing the data. Defaults to \code{sum},
#' which is ideal for count data. For proportioanl data \code{mean} or \code{median}
#' might be more appropriate. Any function that takes a vector of numbers and returns
#' a single numeric value can be used.
#' 
#' @return The default method returns a \code{matrix} object.
#' 
#' The Dataset method returns a \code{Dataset} object when Dat includes a Tax element
#' (see \code{create_dataset}); when the Tax element is missing it returns a
#' \code{matrix} object.
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
#' 
#' @examples 
#' 
pool_samples <- function(x, ...) UseMethod("pool_samples")

#' @rdname pool_samples
#' @method pool_samples default
#' @export
pool_samples.default <- function(x, groups, FUN = sum){
  if(class(x) != "matrix")
    stop("ERROR: Tab must be a matrix object",call.=TRUE)
  if(ncol(x) != length(groups))
    stop("ERROR: Number of columns in Tab and length of groups do not match",call.=TRUE)
  
  groups <- factor(groups)
  res <- collapse_matrix(x = x, groups = groups,
                         dim = 2, FUN = FUN)
  # cat(class(res), "\n")
  
  return(res)
}

#' @rdname pool_samples
#' @method pool_samples Dataset
#' @export
pool_samples.Dataset <- function(x, groups, FUN = sum){
  # Error checking
  if(class(x) != "Dataset")
    stop("ERROR: Dat must be a Dataset object",call.=TRUE)
  if(length(groups) > 1)
    stop("ERROR: groups must be a single value",call.=TRUE)
  if(class(groups) == "character"){
    if(all(groups != colnames(x$Map)))
      stop("ERROR: specified column name not found on Dat$Map",call.=TRUE)
  }else if(class(groups) == "numeric"){
    if(groups > ncol(x$Map))
      stop("ERROR: specified dimension out of range for Dat$Map",call.=TRUE)
  }else{
    stop("ERROR: groups must be either a character or a numeric value indicating which column to use from Dat$Map",call.=TRUE)
  }
  
  res <- pool_samples.default(x = x$Tab, 
                              groups = x$Map[,groups],
                              FUN = FUN)
  
  if(!is.null(x$Tax)){
    res <- create_dataset(Tab = res, Tax = x$Tax)
  }
  
  return(res)
}