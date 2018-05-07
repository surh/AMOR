#' Fit Negative Binomial Generalized Linear Models on a matrix
#' of observations
#' 
#' Takes a matrix of observations, and fits the specified Negative
#' Binomial GLM on each species, returning the matrix of coefficients.
#' 
#' This function will take each species (row) in the abundance matrix
#' and fit the specified (via the formula option) Negative binomial GLM
#' independently on each one, using the glm.nb function from then MASS
#' package.
#' 
#' @param x Either a matrix of abundances with samples as columns and
#' species as rows or a Dataset object
#' @param Map A data.frame containing the variables to be modelled as
#' columns and samples as rows. The rows should be named with sample
#' IDs and must correspond to the column names from x if an abundance
#' matrix was passed
#' @param formula A formula specifying the model to be fit. Only the
#' right hand side of the equation must be passed, and anything on the
#' left side wild be silently ignored
#' @param response.name String indicating the name to be used for the
#' response (dependent) variable in the GLM
#' @param verbose Logical value. If true the taza name is printed while
#' the fit is in progress.
#' @param ... Other parameters for the \link{glm.nb} function like control.
#'
#' @return Returns a matrix.glm object which is a list containing the following elements:
#' \describe{
#'   \item{coefficients}{A matrix of coefficients for the glm fit. The matix has dimensions S x p, where S is the number of species in the abundance matrix, and p the number of parameters in the specified GLM. Each row corresponds to an independent fit of the model with glm() function and the specified formula}
#'   \item{SE}{Similar to Coef, but this matrix contains the estimated Standard Errors (ie. the Standard deviation of the coefficients).}
#'   \item{AIC}{A named vector containing the Akaike Information Criteria (AIC) of each independent fit. The names of the vector correspond to the species IDs in the input abundance matrix}
#'   \item{theta}{A vector of theta (ie. the overdispersion parameter) for each taxon}
#'   \item{SE.theta}{A vector of standard errorrs for the overdispersion parameter.}
#'   \item{call}{Function call via match.call}
#'   \item{family}{GLM model family. An object of class \code{family}}
#'   \item{X}{Design matrix of the model fit.}
#' }
#' @export
#' @author Sur Herrera Paredes
#' @seealso \link{matrix_glm}
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' data(Rhizo.tax)
#' Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)
#' 
#' m1.nb <- matrix_glmNB(x=Dat$Tab,Map=Dat$Map,
#'                       formula=~accession + plate)
#' m1.nb <- matrix_glmNB(Dat,formula=~accession + plate,
#'                       control=glm.control(epsilon=1e-5,maxit=25),
#'                       verbose=TRUE)
matrix_glmNB <- function(x,...) UseMethod("matrix_glmNB")

#' @rdname matrix_glmNB
#' @method matrix_glmNB default
matrix_glmNB.default <- function(x=NULL, Map=NULL, formula=NULL,
                               response.name="Count",
                               verbose=FALSE, ...){

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
  
  Res <- list(coefficients = Beta,SE = SE,VCOV = VCOV, AIC = AIC,call = match.call(),
              family = negative.binomial(theta=theta[ row.names(x) ]),X = X,
              theta = theta,SE.theta = SE.theta)
  #Res <- list(coefficients=Beta,SE=SE,AIC=AIC,family=family)
  class(Res) <- c("matrix.glmNB","matrix.glm")
  return(Res)
}

#' @rdname matrix_glmNB
#' @method matrix_glmNB Dataset
matrix_glmNB.Dataset <- function(x,formula=NULL,response.name="Count",verbose=FALSE,...){
  m1 <- matrix_glmNB.default(x=x$Tab,Map=x$Map,formula=formula,
                           response.name=response.name,verbose=verbose,...)
  m1$call=match.call()
  return(m1)
}
