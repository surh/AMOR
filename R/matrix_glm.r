#' Fit Generalized Linear Models on a matrix of observations
#' 
#' Takes a matrix of observations, and fits the specified GLM
#' on each species, returning the matrix of coefficients.
#' 
#' This function will take each species (row) in the abundance
#' matrix and fit the specified (cia the formula option) GLM
#' independently on each one.
#'
#' @param x Either a matrix of abundances with samples as columns and species as rows or a Dataset object
#' @param Map A data.frame containing the variables to be modelled as columns and samples as rows. The rows should be named with sample IDs and must correspond to the column names from x if an abundance matrix was passed
#' @param formula A formula specifying the model to be fit. Only the right hand side of the equation must be passed, and anything on the left side wild be silently ignored
#' @param family The GLM fmaily to be used. It must be an object of class family. use ?family for more help
#' @param response.name String indicating the name to be used for the response (dependent) variable in the GLM
#' @param verbose Logical value. If true the taza name is printed while the fit is in progress.
#' @param ... Other parameters for the glm() function like control.
#'
#' @return Returns a matrix.glm object which is a list containing the following elements:
#' \describe{
#'   \item{coefficients}{A matrix of coefficients for the glm fit. The matix has dimensions S x p, where S is the number of species in the abundance matrix, and p the number of parameters in the specified GLM. Each row corresponds to an independent fit of the model with glm() function and the specified formula}
#'   \item{SE}{Similar to Coef, but this matrix contains the estimated Standard Errors (ie. the Standard deviation of the coefficients).}
#'   \item{AIC}{A named vector containing the Akaike Information Criteria (AIC) of each independent fit. The names of the vector correspond to the species IDs in the input abundance matrix}
#'   \item{call}{Function call via match.call}
#'   \item{family}{GLM model family. An object of class \code{family}}
#'   \item{X}{Design matrix of the model fit.}
#' }
#' @author Sur Herrera Paredes
#' @export
matrix_glm <- function(x, Map, formula, family, response.name,
                       verbose, ...) UseMethod("matrix_glm")

#' @rdname matrix_glm
#' @method matrix_glm default
matrix_glm.default <- function(x=NULL,Map=NULL,
                               formula=NULL,family=poisson(link="log"),
                               response.name="Count",verbose=FALSE, ...){
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
  VCOV <- rep(list(matrix(NA)), times = length(row.names(x)))
  names(VCOV) <- row.names(x)
  op <- options(warn=1)
  for(taxa in row.names(x)){
    #taxa <- "79"
    #taxa <- "OTU_14834"
    
    count <- x[taxa,]
    Map[,response.name] <- as.numeric(count)
    if(verbose) cat(taxa,"\n") 
    m1 <- list(coefficients=NA)
    m1 <- tryCatch(glm(f1,data=Map,family=family,...),
             error=function(e){list(coefficients=NA)})
    
    if(!is.na(m1$coefficients[1])){
      # Get fixed effect coefficients
      Beta[taxa,] <- m1$coefficients
      VCOV[[taxa]] <- vcov(m1)
      SE[taxa,] <- sqrt(diag(VCOV[[taxa]]))[names(m1$coefficients)]
      AIC[taxa] <- AIC(m1)
    }

    #cat("...Done\n")
  }
  options(op)
  
  Res <- list(coefficients = Beta,SE = SE,VCOV = VCOV,
              AIC = AIC,call = match.call(),family = family,X = X)
  #Res <- list(coefficients=Beta,SE=SE,AIC=AIC,family=family)
  class(Res) <- "matrix.glm"
  return(Res)
}

#' @rdname matrix_glm
#' @method matrix_glm Dataset
matrix_glm.Dataset <- function(x,formula=NULL,
                               family=poisson(link="log"),
                               response.name="Count",
                               verbose=FALSE,...){
  m1 <- matrix_glm.default(x=x$Tab, Map=x$Map, formula=formula,
                           family=family,
                           response.name=response.name,
                           verbose=verbose, ...)
  m1$call=match.call()
  return(m1)
}

#' Summarize matrix GLM results
#'
#'Takes a matrix.glm or bootglm object and reformats the results
#'
#' @param object A matrix.glm or bootglm object.
#' @param sortby A string indicating how to sort the coefficients.
#' Either by "Vatiable" or "Taxon"
#' 
#' @return Returns a summary.matrix.glm or summary.bootglm
#' object which is a list containing the following elements:
#' \describe{
#' \item{coefficients}{Coefficient table. For matrix.glm objects it
#' includes standard errors, z-values and p-values. NOTE: Since only
#' p-values from z-tests are provided at this point, they are only
#' appropriate for binomial, poisson and negative.binomial models.
#' For boot.glm objects in includes the 2.5% and 97.5% percentiles
#' and the bootstrap confidence for each parameter (Pr(param > 0)).}
#' \item{call}{Original matrix_glm() or bootstrap_glm() call}}
#' @author Sur Herrera Paredes
#' @export
summary.matrix.glm <- function(object,sortby="Variable"){
  
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

#' @export
print.summary.matrix.glm<-function(x,...){
  cat("Call:\n")
  print(x$call)
  
  row.names(x$coefficients) <- paste(as.character(x$coefficients$Taxon),
                                     as.character(x$coefficients$Variable),
                                     sep="__")
  x$coefficients$Taxon <- NULL
  x$coefficients$Variable <- NULL
  
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE, digits=3)
  
}

