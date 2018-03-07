collapse_matrix <- function(x,groups,dim=1,FUN=sum){
  #   x <- as.matrix(Rhizo)
  #   groups <- rep(c(1,2),length.out=nrow(x))
  #   dim <- 1
  #   groups <- rep(c(1,2),length.out=ncol(x))
  #   dim <- 2
  #   FUN <- sum
  
  if(class(x) != "matrix")
    stop("ERROR: Tab must be a matrix ",call.=TRUE)
  if(dim != 1 & dim != 2 )# NOTE: Migth eventually allow for arbitrary dimensiton arrays.
    stop("ERROR: dim must be 1 or 2")
  if(length(groups) != length(dimnames(x)[[dim]])){
    stop("ERROR: length of groups does not correspond to number of vectors",call.=TRUE)
  }
  
  FUN <- match.fun(FUN)
  groups <- factor(groups)  
  Tab.collapsed <- apply(x,((dim %% 2) + 1),tapply,
                         INDEX = groups, FUN)
  if(dim == 2){
    Tab.collapsed <- t(Tab.collapsed)
  }
  return(Tab.collapsed)
}