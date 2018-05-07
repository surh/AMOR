#' Principal Coordinates Analysis.
#' 
#' Function that performs Principal Coordinates Analysis on abundance matrix.
#' 
#' @param x Distance matrix, must be a \code{dist} object, or a \code{Dataset}
#' object. See \code{\link{create_dataset}} for more information.
#' @param dim Number of dimensions to return.
#' @param distfun Function that calculates distance matrix for abundance
#' table in \code{Dataset} object.
#' 
#' @details This is the same as function \code{\link{pco}} from the labdsv,
#' but it includes a \code{Dataset} method.
#' 
#' @return A PCO object as defined by the \code{\link{pco}} function.
#' 
#' @author Sur from Dangl Lab.
#' 
#' @seealso \code{\link{create_dataset}}, \code{\link{pca}}, \code{\link{PCA}},
#' \code{\link{pco}}, \code{\link{dsvdis}}
#' 
#' @export
#' 
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' Dat <- create_dataset(Rhizo,Rhizo.map)
#' 
#' # distfun <- function(x) vegan::vegdist(x,method="bray") #requires vegan package
#' distfun <- dist
#' 
#' Dat.pco <- PCO(Dat,dim=2,distfun=distfun)
#' summary(Dat.pco)
#' plotgg(Dat.pco)
#' plotgg(Dat.pco,shape="fraction",point_size=3)
#' plotgg(Dat.pco,shape="fraction",col="accession",point_size=4)
PCO <- function(x, dim, distfun) UseMethod("PCO")

#' @rdname PCO
#' @method PCO default
#' @export
PCO.default <- function(x, dim=3){
  # Taken from labdsv
  res <- cmdscale(x, k = dim, eig = TRUE)
  res$Map <- NULL
  res$Tax <- NULL
  colnames(res$points) <- paste("PCo",1:dim,sep = "")
  class(res) <- "PCO"
  return(res)
}

#' @rdname PCO
#' @method PCO Dataset
#' @export
PCO.Dataset <- function(x, dim=3, distfun=dist){
  mat <- x$Tab
  mat <- t(mat)
  mat.dist <- distfun(mat)
  mat.pco <- PCO.default(mat.dist,dim=dim)
  mat.pco$Map <- x$Map
  mat.pco$Tax <- x$Tax
  return(mat.pco)
}

#' @export
summary.PCO <- function(object){
  #object <- Dat.pco
  ncomponents <- ncol(object$points)
  components <- paste("PCo",1:ncomponents,sep = "")
  percvar <- 100*(object$eig[ object$eig >= 0 ]) / sum(object$eig[ object$eig >= 0 ])
  cumvar <- cumsum(percvar)
  
  vartab <- data.frame(Component = components,
                       Var.explained = percvar[1:ncomponents],
                       Cumulative = cumvar[1:ncomponents])
  
  sum.pco <- list(vartab = vartab,
                  ncomponents = ncomponents)
  class(sum.pco) <- "summary.PCO"
  
  return(sum.pco)
}

#' @export
print.summary.PCO <- function(x,digits = 2, n = 5){
  cat("Principal Coordinate Analysis:\n")
  cat("\t",x$ncomponents, " Components\n\n")
  
  tab <- x$vartab[ 1:min(n,nrow(x$vartab)), ]
  tab[ ,2:3 ] <- round(tab[,2:3],digits = digits)
  print(tab)
}