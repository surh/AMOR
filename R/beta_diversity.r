#' Beta diversity
#' 
#' Calculate beta diversity from abundance matrix
#'
#' @param X Matrix with of abundances
#' @param method Distance method to be used. One of: "Bray-Curtis"
#' @param sample_dim Dimension that contains samples, 1 for rows
#' and 2 for columns. The other dimension should contain the taxa
#'
#' @return A \code{dist} object
#' 
#' @seealso \link{dist}
#' 
#' @export
#' 
#' @author Sur Herrera Paredes
#'
#' @examples
#' mat <- matrix(1:10,ncol = 2)
#' beta_diversity(mat)
beta_diversity <- function(X,method="Bray-Curtis",sample_dim=2){
  X <- as.matrix(X)
  
  if(method == "Bray-Curtis"){
    if(sample_dim == 1){
      X.t <- t(X)
    }else if(sample_dim == 2){
      X.t <- X
    }else{
      stop("Invalid sample_dim",call=TRUE)
    }
    X.d <- apply(X,sample_dim,function(vec,X.t){  
      num <- colSums(abs(vec - X.t))
      den <- colSums(vec + X.t)
      d <- num/den  
      d
    },X.t=X.t)
    
    X.d <- as.dist(X.d)
  }else{
    stop("Unknown method",call=TRUE)
  }
  
  return(X.d)
}
