matrix_glmNB <- function(x,...) UseMethod("matrix_glmNB")

matrix_glmNB.default <- function(x=NULL,Map=NULL,formula=NULL,
                               response.name="Count",verbose=FALSE,...){
  # Test data
#     formula <- formula(~ frac + gen + depthK)
#     x <- Dat$Tab
#     Map <- Dat$Map
#     response.name <- "Count"
  
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
  theta <- rep(NA,times=length(row.names(x)))
  SE.theta <- rep(NA,times=length(row.names(x)))
  names(AIC) <- row.names(x)
  names(theta) <- row.names(x)
  names(SE.theta) <- row.names(x)
  VCOV <- rep(list(matrix(NA)), times = length(row.names(x)))
  names(VCOV) <- row.names(x)
  op <- options(warn=1)
  for(taxa in row.names(x)){
    #taxa <- "79"
    count <- x[taxa,]
    Map[,response.name] <- as.numeric(count)
    if(verbose) cat(taxa,"\n") 
    m1 <- list(coefficients=NA)
    m1 <- tryCatch(glm.nb(f1,data=Map,...),
                   error=function(e){list(coefficients=NA)})
    
    if(!is.na(m1$coefficients[1])){
      # Get fixed effect coefficients
      Beta[taxa,] <- m1$coefficients
      VCOV[[taxa]] <- vcov(m1)
      SE[taxa,] <- sqrt(diag(VCOV[[taxa]]))[names(m1$coefficients)]
      AIC[taxa] <- AIC(m1)
      theta[taxa] <- m1$theta
      SE.theta[taxa] <- m1$SE.theta
    }
    
    #cat("...Done\n")
  }
  options(op)
  
  Res <- list(coefficients=Beta,SE=SE,VCOV = VCOV, AIC=AIC,call=match.call(),
              family=negative.binomial(theta=theta[ row.names(x) ]),X=X,
              theta=theta,SE.theta=SE.theta)
  #Res <- list(coefficients=Beta,SE=SE,AIC=AIC,family=family)
  class(Res) <- c("matrix.glmNB","matrix.glm")
  return(Res)
}

matrix_glmNB.Dataset <- function(x,formula=NULL,response.name="Count",verbose=FALSE,...){
  m1 <- matrix_glmNB.default(x=x$Tab,Map=x$Map,formula=formula,
                           response.name=response.name,verbose=verbose,...)
  m1$call=match.call()
  return(m1)
}
