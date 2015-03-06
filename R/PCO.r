PCO <- function(...) UseMethod("PCO")

PCO.default <- function(x,dim=3){
  # Taken from labdsv
  res <- cmdscale(x, k = dim, eig = TRUE)
  res$Map <- NULL
  res$Tax <- NULL
  class(res) <- "PCO"
  return(res)
}

PCO.Dataset <- function(Dat,dim=3,distfun=dist){
  mat <- Dat$Tab
  mat <- t(mat)
  mat.dist <- distfun(mat)
  mat.pco <- PCO.default(mat.dist,dim=dim)
  mat.pco$Map <- Dat$Map
  mat.pco$Tax <- Dat$Tax
  return(mat.pco)
}
