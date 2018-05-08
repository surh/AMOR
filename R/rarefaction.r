#' Rarefaction
#' 
#' Performs rarefaction on count table.
#' 
#' Modified from the rrarefy function of vegan-2.0-4.
#' 
#' @param x Either a matrix object or a Dataset object.
#' If a matrix is passed, samples must be the columns
#' and rows the species.
#' @param sample Sample size, either a single integer
#' value or a vector of the same length the number of
#' columns (samples) of x.
#'
#' @return The defualt method returns a matrix with
#' the rarefied counts. The Dataset method returns a
#' Dataset object where the Tab element is the rarefied
#' counts.
#' @export
#' @author Sur Herrera Paredes
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' Dat <- create_dataset(Rhizo,Rhizo.map)
#' set.seed(712) # Always specify and save the seed before rarefaction
#' rarefaction(x=Dat$Tab,sample=100)
#' set.seed(712)
#' rarefaction(x=Dat,sample=100)$Tab
rarefaction <- function(x, sample) UseMethod("rarefaction")

#' @rdname rarefaction
#' @method rarefaction default
#' @export
rarefaction.default <- function (x, sample){
  if (!identical(all.equal(x, round(x)), TRUE)){
    stop("function is meaningful only for integers (counts)")
  }
  if(class(x) != "matrix"){
    stop("ERROR: function only works on matrix objects",call.=TRUE)
  }
  if (length(sample) > 1 && length(sample) != ncol(x)){
    stop(gettextf("length of 'sample' and number of columns of 'x' do not match"))
  }
  sample <- rep(sample, length = ncol(x))
  rownames(x) <- rownames(x, do.NULL = FALSE)
  nm <- rownames(x)
  
  for (i in 1:ncol(x)) {
    col <- sample(rep(nm, times = x[,i]), sample[i])
    col <- table(col)
    ind <- names(col)
    x[,i] <- 0
    x[ind,i] <- col
  }
  return(x)
}

#' @rdname rarefaction
#' @method rarefaction Dataset
#' @export
rarefaction.Dataset <- function (x, sample){
  Tab <- rarefaction.default(x=x$Tab,sample=sample)
  x <- create_dataset(Tab=Tab,Map=x$Map,Tax=x$Tax)
  return(x)
}