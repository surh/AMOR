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

summary.PCO <- function(object){
  object <- Dat.pco
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

print.summary.PCO <- function(x,digits = 2, n = 5){
  cat("Principal Coordinate Analysis:\n")
  cat("\t",x$ncomponents, " Components\n\n")
  
  tab <- x$vartab[ 1:min(n,nrow(x$vartab)), ]
  tab[ ,2:3 ] <- round(tab[,2:3],digits = digits)
  print(tab)
}


