rarefaction <- function(x,sample) UseMethod("rarefaction")

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

rarefaction.Dataset <- function (x, sample){
  Tab <- rarefaction.default(x=x$Tab,sample=sample)
  x <- create_dataset(Tab=Tab,Map=x$Map,Tax=x$Tax)
  return(x)
}

