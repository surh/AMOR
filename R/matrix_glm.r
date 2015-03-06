matrix_glm <- function(x,...) UseMethod("matrix_glm")

matrix_glm.default <- function(x=NULL,Map=NULL,formula=NULL,family=poisson(link="log"),
                   response.name="Count",verbose=FALSE,...){
  # Test data
#   formula <- formula(~ frac + gen + depthK)
#   x <- Dat$Tab
#   Map <- Dat$Map
#   control <- list(maxit=100)
#   family <- poisson(link="log")
#   response.name <- "Count"
  
#     formula <- formula(~ frac + gen + depthK)
#     x <- Tab.bs
#     Map <- Dat.bs
#     Map$frac <- factor(Map$frac,levels=c("R","E"))
#     control <- list(maxit=100)
#     family <- poisson(link="log")
#     response.name <- "Count"
    
  # Check names
  if(any(colnames(x) != row.names(Map))){
    stop("ERROR: sample names in abundance table don't match sample names in mapping file",call.=TRUE)
  }
  if(response.name %in% names(Map)){
    stop("ERROR: response name already exists in mapping file",call.=TRUE)
  }
  
  X <- model.matrix(formula,data=Map)
  coef_names <- colnames(X)
  f1 <- update(formula,paste(response.name," ~ ."))
  
  # Fit
  Beta <- matrix(NA,ncol=length(coef_names),nrow=length(row.names(x)),
                 dimnames=list(row.names(x),coef_names))
  SE <- Beta
  AIC <- rep(NA,times=length(row.names(x)))
  names(AIC) <- row.names(x)
  op <- options(warn=1)
  for(taxa in row.names(x)){
    #taxa <- "79"
    count <- x[taxa,]
    Map[,response.name] <- as.numeric(count)
    if(verbose) cat(taxa,"\n") 
    m1 <- list(coefficients=NA)
    m1 <- tryCatch(glm(f1,data=Map,family=family,...),
             error=function(e){list(coefficients=NA)})
    
    if(!is.na(m1$coefficients[1])){
      # Get fixed effect coefficients
      Beta[taxa,] <- m1$coefficients
      SE[taxa,] <- sqrt(diag(vcov(m1)))[names(m1$coefficients)]
      AIC[taxa] <- AIC(m1)
    }

    #cat("...Done\n")
  }
  options(op)
  
  Res <- list(coefficients=Beta,SE=SE,AIC=AIC,call=match.call(),family=family,X=X)
  #Res <- list(coefficients=Beta,SE=SE,AIC=AIC,family=family)
  class(Res) <- "matrix.glm"
  return(Res)
}

matrix_glm.Dataset <- function(x,formula=NULL,family=poisson(link="log"),response.name="Count",verbose=FALSE,...){
  m1 <- matrix_glm.default(x=x$Tab,Map=x$Map,formula=formula,family=family,
                           response.name=response.name,verbose=verbose,...)
  m1$call=match.call()
  return(m1)
}

summary.matrix.glm <- function(object,sortby="Variable",...){
  # Test data
#   object <- m1.bs$orig
#   sortby <- "Variable"
  
  if(sortby != "Taxon" && sortby != "Variable"){
    stop("ERROR: You can only sort by Taxon or Variable",call.=TRUE)
  }
  
  TAB <- melt(object$coefficients,value.name="Estimate",varnames=c("Taxon","Variable"))
  se <- melt(object$SE,value.name="StdErr",varnames=c("Taxon","Variable"))[,"StdErr"]
  TAB$StdErr <- se
  
  TAB$z.value <- TAB$Estimate / TAB$StdErr
  TAB$p.value <- 2*pnorm(-abs(TAB$z.value))
  
  if(sortby == "Taxon"){
    TAB <- TAB[order(TAB[,sortby]),]
  }
  
  Res <- list(coefficients=TAB,call=object$call)
  class(Res) <- "summary.matrix.glm"
  return(Res)
}

print.summary.matrix.glm<-function(x,...){
  # Test
#   x <- Res
  cat("Call:\n")
  print(x$call)
  
  row.names(x$coefficients) <- paste(as.character(x$coefficients$Taxon),as.character(x$coefficients$Variable),sep="__")
  x$coefficients$Taxon <- NULL
  x$coefficients$Variable <- NULL
  
  printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE,digits=3)
  
}

