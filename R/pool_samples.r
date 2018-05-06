#' Pool samples
#' 
#' Pools abundances from samples according to a grouping factor.
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
#' @param groups For the default and Dataset methods this can be either a vector
#' or factor specifying to which group each sample belongs. Vectors will be
#' converted to a factor with the \code{factor} function. It must be of the same length
#' as \code{ncol(x)}. For the Dataset method, it can also be a single character string or
#' numeric value indicating the column from x$Map to be used as a grouping factor.
#' @param FUN Function to apply when collapsing the data. Defaults to \code{sum},
#' which is ideal for count data. For proportioanl data \code{mean} or \code{median}
#' might be more appropriate. Any function that takes a vector of numbers and returns
#' a single numeric value can be used.
#' @param return.dataset Logical, if TRUE returns a dataset
#' 
#' @return The default method returns a \code{matrix} object, unless return.dataset
#' is TRUE, in which case it returns a Dataset.
#' 
#' The Dataset method returns a \code{Dataset} object when Dat includes a Tax element
#' (see \code{create_dataset}); when the Tax element is missing it returns a
#' \code{matrix} object. If return.dataset is TRUE, return a dataset
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
#' 
#' @examples 
#' library(AMOR)
#' data(Rhizo)
#' data(Rhizo.map)
#' data(Rhizo.tax)
#' Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)
#' 
#' # The following returns a numeric matrix
#' Collapsed1 <- pool_samples(x = Dat$Tab,groups = Dat$Map$fraction)
#' 
#' # The following returns a Dataset
#' Collapsed2 <- pool_samples(x = Dat,groups = "fraction")
#' 
#' # You can also directly pass a grouping factor to the Dataset method
#' Collapsed3 <- pool_samples(x = Dat,groups = Dat$Map$fraction)
#' 
#' # A way to calculate the overall counts per taxa
#' res <- pool_samples(Dat$Tab, groups = rep("all", length.out = ncol(Rhizo)))
pool_samples <- function(x, ...) UseMethod("pool_samples")

#' @rdname pool_samples
#' @method pool_samples default
#' @export
pool_samples.default <- function(x, groups, FUN = sum, return.dataset = FALSE){
  if(class(x) != "matrix")
    stop("ERROR: Tab must be a matrix object",call.=TRUE)
  if(ncol(x) != length(groups))
    stop("ERROR: Number of columns in Tab and length of groups do not match",call.=TRUE)
  
  groups <- factor(groups)
  res <- collapse_matrix(x = x, groups = groups,
                         dim = 2, FUN = FUN)
  # cat(class(res), "\n")
  if (return.dataset)
    res <- create_dataset(res)
  
  return(res)
}

#' @rdname pool_samples
#' @method pool_samples Dataset
pool_samples.Dataset <- function(x, groups, FUN = sum, return.dataset = FALSE){
  # Check object class
  if(class(x) != "Dataset")
    stop("ERROR: Dat must be a Dataset object",call.=TRUE)
  
  # Check grouping factor
  if(length(groups) == length(samples(x))){
    groups <- factor(groups)
  }else if(length(groups) == 1){
    if(class(groups) == "character"){
      if(all(groups != colnames(x$Map)))
        stop("ERROR: specified column name not found on Dat$Map",call.=TRUE)
    }else if(class(groups) == "numeric"){
      if(groups > ncol(x$Map))
        stop("ERROR: specified dimension out of range for Dat$Map",call.=TRUE)
    }else{
      stop("ERROR: groups must be either a character or a numeric value indicating which column to use from Dat$Map",call.=TRUE)
    }
    
    groups <- factor(x$Map[ , groups ])
  }else{
    stop("ERROR: groups must be either a vector or factor of length length(samples(x)) indicating how to group the samples, or a single numeric or character value indicating a variable in x$Map to use as grouping factor",call.=TRUE)
  }
    
  res <- pool_samples.default(x = x$Tab, 
                              groups = groups,
                              FUN = FUN)
  if(return.dataset){
    res <- create_dataset(Tab = res, Tax = x$Tax)
  }
  if(!is.null(x$Tax)){
    res <- create_dataset(Tab = res, Tax = x$Tax)
  }
  
  return(res)
}