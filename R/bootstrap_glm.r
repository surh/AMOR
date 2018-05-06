#' Bootstrap estimation of Generalized Linear Models on a matrix of observations
#' 
#' Uses bootstrap to estimate variance of \link{matrix_glm}
#' estimates. Takes a matrix of observations, and fits the specified
#' GLM on each species. It perfomrs N bootstrap pesudoreplicates and
#' returns the matrix of coefficients, for the original dataset as well
#' as bootstrap estimates. It uses glm_matrix() to fit each of the pseudoreplicates.
#' 
#' This works by calling matrix_glm() on the original data, and each pf the N
#' bootstrap pseudoreplicates generated.
#'
#' @param x Either an abundance matrix with samples as columns
#' or a Dataset object
#' @param Map A mapping file with sample metadata if x is a matrix. A data.frame
#' containing the variables to be modelled as columns and samples as rows.
#' The rows should be named with sample IDs and must correspond to the column names
#' from x if an abundance matrix was passed.
#' @param formula The formula for the GLM model.Only the right hand side of the equation
#' must be passed, and anything on the left side wild be silently ignored.
#' @param family The model family. See help of \link{glm} and \link{family}
#' for details
#' @param response.name The name to give to the response variable
#' @param N The number of bootstrap pseudoreplicates
#' @param verbose Logical. If true each call to matrix_glm will print
#' progress information.
#' @param ... Extra parameters for \link{matrix_glm}, which will in turn
#' pass them to \link{glm}. Typically it is seful to specify the 'control'
#' option.
#'
#' @return A bootglm object which is a list containing the following elements:
#' \describe{
#' \item{orig}{A \code{matrix.glm} object containing the output of \code{matrix.glm()}
#' on the original data.}
#' \item{boot}{A list containing elements: 1) coefficients, a three-dimensioanl
#' array containing the matrix of coefficients for each bootstrap pseudoreplicate;
#' and 2) SE, a three-dimensional array containing the matrix of standard errors for
#' each bootstrap pseudoreplicate.}
#' \item{coefficients}{A matrix of coefficients calculated as the mean accross each
#' bootstrap pseudoreplicates.}
#' \item{call}{Function call via match.call().}
#' \item{family}{GLM model family. An object of class \code{family}.}
#' \item{pseudoreplicates}{Matrix of dimensions (number iof samples) x (number of
#' pseudoreplicates), where the i-column gives the sample ID, of the samples used for
#' the i-th pseudoreplicate.}
#' \item{N}{Number of bootstrap pseudoreplicates performed.}
#' }
#' @author Sur Herrera Paredes
#' @export
bootstrap_glm <- function(x,...) UseMethod("bootstrap_glm")

#' @rdname bootstrap_glm
#' @method bootstrap_glm default
bootstrap_glm.default <- function(x,Map,formula,family=poisson(link="log"),
                                  response.name="Count",N=100,verbose=FALSE,...){
  
  cat("Fitting model with orignal data...\n")
  m1 <- matrix_glm(x=x, Map=Map, formula=formula,
                   family=family,
                   response.name=response.name,
                   verbose=verbose, ...)
  
  res.template <- array(dim=c(dim(m1$coefficients),N))
  row.names(res.template) <- row.names(m1$coefficients)
  colnames(res.template) <- colnames(m1$coefficients)
  
  Res <- list(orig=m1,
              boot=list(coefficients=res.template,
                        SE=res.template))
  res.template <- NULL
    
  cat("Fitting bootstrap pseudoreplicates...\n")
  SAMPLES <- matrix(NA,nrow=dim(x)[2],ncol=N)
  for(i in 1:N){
    #i <- 1
    if(i %% 50 == 0) cat("Replicate ",i,"\n")
    pseudorep <- sample(colnames(x),replace=T)
    SAMPLES[,i] <- pseudorep
    Tab.bs <- x[,pseudorep]
    Dat.bs <- Map[colnames(Tab.bs),]
    Dat.bs <- droplevels(Dat.bs)
    colnames(Tab.bs) <- row.names(Dat.bs)
    
    # Remove categorical variables with only one level
    minus_levs <- colnames(Dat.bs)[lapply(lapply(Dat.bs,levels),length) == 1]
    form.char <- as.character(formula)
    right.side <- form.char[length(form.char)]
    included.vars <- strsplit(right.side,split="[+]")[[1]]
    included.vars <- gsub(pattern=" ",replacement="",x=included.vars)
    final.vars <- included.vars[!(included.vars %in% minus_levs)]
    formula <- paste(final.vars,collapse=" + ")
    formula <- paste("~",formula)
    formula <- formula(formula)
    
    m1 <- matrix_glm(x=Tab.bs, Map=Dat.bs,
                     formula=formula, family=family,
                     response.name=response.name,
                     verbose=verbose,c...)
    
    Res$boot$coefficients[,colnames(m1$coefficients),i] <- m1$coefficients
    Res$boot$SE[,colnames(m1$SE),i] <- m1$SE
  }
  
  Res$call <- match.call()
  Res$family <- family
  Res$pseudoreplicates <- SAMPLES
  Res$N <- N
  Res$coefficients <- apply(Res$boot$coefficients,c(1,2), mean, na.rm=TRUE)
  class(Res) <- "bootglm"
  
  return(Res)
}

#' @rdname bootstrap_glm
#' @method bootstrap_glm Dataset
bootstrap_glm.Dataset <- function(x,formula, family=poisson(link="log"),
                                  response.name="Count",
                                  N=100, verbose = FALSE,...){
  m1 <- bootstrap_glm.default(x=x$Tab, Map=x$Map, 
                              formula=formula,
                              family=family,
                              response.name=response.name,
                              N=N, ...)
  m1$call <- match.call()
  return(m1)
}


#' Summarizing bootstrap_glm model fits
#'
#' @param object A bootglm object
#' @param sortby How to sort the summarized results
#'
#' @return A summary.bootglm object
#' @author Sur Herrera Paredes
#' @export
summary.bootglm <- function(object,sortby="Variable",...){
  if(sortby != "Taxon" && sortby != "Variable"){
    stop("ERROR: You can only sort by Taxon or Variable",call.=TRUE)
  }
  
  TAB <- melt(object$coefficients, value.name="Estimate", varnames=c("Taxon","Variable"))
  quants <- apply(object$boot$coefficients,c(1,2),
                  function(vector){
                    x <- quantile(vector, probs=c(0.025,0.975), na.rm=TRUE) 
                    n <- sum(!is.na(vector))
                    conf <- sum(vector > 0,na.rm=TRUE) / n
                    x <- c(x,N=n,conf=conf)
                    x
                  })
  
  TAB[,"2.5%"] <- melt(quants[1,,],value.name="2.5%",varnames=c("Taxon","Variable"))[,3]
  TAB[,"97.5%"] <- melt(quants[2,,],value.name="2.5%",varnames=c("Taxon","Variable"))[,3]
  TAB[,"N"] <- melt(quants[3,,],value.name="2.5%",varnames=c("Taxon","Variable"))[,3]
  TAB[,"Pr(>0)"] <- melt(quants[4,,],value.name="2.5%",varnames=c("Taxon","Variable"))[,3]
  
  if(sortby == "Taxon"){
    TAB <- TAB[order(TAB[,sortby]),]
  }
  
  Res <- list(coefficients=TAB,call=object$call)
  class(Res) <- "summary.bootglm"
  return(Res)
}

#' Print summary.bootglm object
#'
#' @param x A summary.bootglm object
#' 
#' @author Sur Herrera Paredes
#' @export
print.summary.bootglm<-function(x){
  # Test
  #   x <- Res
  cat("Call:\n")
  print(x$call)
  
  row.names(x$coefficients) <- paste(as.character(x$coefficients$Taxon),
                                     as.character(x$coefficients$Variable),sep="__")
  x$coefficients$Taxon <- NULL
  x$coefficients$Variable <- NULL
  x$N <- format(x$N,digits=2)
  
  print(x$coefficients)
}