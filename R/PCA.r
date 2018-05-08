#' Principal Component Analysis
#' 
#' Function that performs principal component analysis on an abundance matrix.
#' 
#' @param x Numeric matrix where samples are columns and rows are species, or a
#' Dataset object, see \code{\link{create_dataset}}.
#' @param cor logical value indicating whether the correlation matrix should
#' be used instead of the covariance matrix.
#' @param dim Number of dimensions to return.
#' 
#' @details This function is the same as function \code{\link{pca}} from the
#' \code{\link{labdsv}} package, but includes a methdod for Dataset objects.
#' 
#' @return A \code{PCA} object. Includes the same attributes as a
#' \code{\link{pca}} object from the \code{\link{labdsv}} package.
#' When the Dataset method is used, it includes two additional slots:
#' \itemize{
#' \item{"Map"}{The Mapping file for the samples.}
#' \item{"Tax"}{The Taxonomic information of the taxa.}}
#' 
#' @author Sur from Dangl Lab.
#' 
#' @seealso  \link{create_dataset}, \link{pca}, \link{PCO}, \link{pco},
#' \link{plotgg.pca}
#' 
#' @export
#' 
#' @examples 
#' data(Rhizo)
#' data(Rhizo.map)
#' Dat <- create_dataset(Rhizo,Rhizo.map)
#' Dat.pca <- PCA(Dat)
#' plotgg(Dat.pca,col="accession",shape="fraction",point_size=4,biplot=TRUE)
#' summary(Dat.pca)
PCA <- function(x, cor, dim) UseMethod("PCA")

#' @rdname PCA
#' @method PCA default
#' @export
PCA.default <- function(x, cor=FALSE,
                        dim=min(nrow(x), ncol(x))){
  # Taken from labdsv
  x <- t(x)
  temp <- prcomp(x,retx=TRUE,center=TRUE,scale=cor)
  res <- list(scores = temp$x[,1:dim],
              loadings = temp$rotation[,1:dim],
              sdev = temp$sdev[1:dim],
              totvar = sum(temp$sdev^2),
              Map = NULL,
              Tax = NULL)
  
  class(res) <- c("PCA")
  return(res)
}

#' @rdname PCA
#' @method PCA Dataset
#' @export
PCA.Dataset <- function(x, cor = FALSE,
                        dim = min(nrow(x$Tab), ncol(x$Tab))){
  mat <- x$Tab
  res <- PCA.default(mat, cor = cor, dim = dim)
  res$Map <- x$Map
  res$Tax <- x$Tax
  class(res) <- c("PCA")
  return(res)
}

#' @export
summary.PCA <- function(object){
  ncomponents <- length(object$sdev)
  components <- paste("PC",1:ncomponents,sep = "")
  percvar <- 100*(object$sdev^2) / object$totvar
  cumvar <- cumsum(percvar)
  
  vartab <- data.frame(Component = components,
                       Var.explained = percvar,
                       Cumulative = cumvar)
  
  sum.pca <- list(vartab = vartab,
                  ncomponents = ncomponents)
  class(sum.pca) <- "summary.PCA"
  
  return(sum.pca)
}

#' @export
print.summary.PCA <- function(x,digits = 2, n = 5){
  cat("Principal Component Analysis:\n")
  cat("\t",x$ncomponents, " Components\n\n")
  
  tab <- x$vartab[ 1:n, ]
  tab[ ,2:3 ] <- round(tab[,2:3],digits = digits)
  print(tab)
}

