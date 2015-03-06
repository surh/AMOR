pool_samples <- function(...) UseMethod("pool_samples")

pool_samples.default <- function(Tab,groups,FUN = sum){
  if(class(Tab) != "matrix")
    stop("ERROR: Tab must be a matrix object",call.=TRUE)
  if(ncol(Tab) != length(groups))
    stop("ERROR: Number of columns in Tab and length of groups do not match",call.=TRUE)
  
  groups <- factor(groups)
  Tab.collapsed <- collapse_matrix(x=Tab,groups=groups,dim=2,FUN=FUN)
  res <- create_dataset(Tab = Tab.collapsed)
  return(res)
}



pool_samples.Dataset <- function(Dat,groups,FUN = sum){
#   groups <- "fraction"
#   FUN <- sum
  
  # Error checking
  if(class(Dat) != "Dataset")
    stop("ERROR: Dat must be a Dataset object",call.=TRUE)
  if(length(groups) > 1)
    stop("ERROR: groups must be a single value",call.=TRUE)
  if(class(groups) == "character"){
    if(all(groups != colnames(Dat$Map)))
      stop("ERROR: specified column name not found on Dat$Map",call.=TRUE)
  }else if(class(groups) == "numeric"){
    if(groups > ncol(Dat$Map))
      stop("ERROR: specified dimension out of range for Dat$Map",call.=TRUE)
  }else{
    stop("ERROR: groups must be either a character or a numeric value indicating which column to use from Dat$Map",call.=TRUE)
  }
  
  res <- pool_samples.default(Tab = Dat$Tab,groups = Dat$Map[,groups],FUN = FUN)
  res <- create_dataset(Tab=res$Tab,Tax=Dat$Tax)
  
  return(res)
}