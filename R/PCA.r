PCA <- function(...) UseMethod("PCA")

PCA.default <- function(x,cor=FALSE,dim=min(nrow(x),ncol(x))){
  # Taken from labdsv
#   res <- labdsv::pca(...)
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

PCA.Dataset <- function(Dat,...){
  mat <- Dat$Tab
  res <- PCA.default(mat,...)
  res$Map <- Dat$Map
  res$Tax <- Dat$Tax
  class(res) <- c("PCA")
  return(res)
}

