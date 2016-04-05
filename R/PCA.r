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

PCA.Dataset <- function(Dat,cor = FALSE, dim = min(nrow(Dat$Tab),ncol(Dat$Tab))){
  mat <- Dat$Tab
  res <- PCA.default(mat, cor = cor, dim = dim)
  res$Map <- Dat$Map
  res$Tax <- Dat$Tax
  class(res) <- c("PCA")
  return(res)
}


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

print.summary.PCA <- function(x,digits = 2, n = 5){
  cat("Principal Component Analysis:\n")
  cat("\t",x$ncomponents, " Components\n\n")
  
  tab <- x$vartab[ 1:n, ]
  tab[ ,2:3 ] <- round(tab[,2:3],digits = digits)
  print(tab)
}

