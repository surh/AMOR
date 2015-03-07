bootstrap_glm <- function(x,...) UseMethod("bootstrap_glm")

bootstrap_glm.default <- function(x,Map,formula,family=poisson(link="log"),
                                  response.name="Count",N=100,verbose=FALSE,...){
  # Test data
#   formula <- formula(~ frac + gen + depthK)
#   x <- Dat$Tab
#   Map <- Dat$Map
#   control <- list(maxit=100)
#   family <- poisson(link="log")
#   response.name <- "Count"
#   N <- 10
#   verbose <- FALSE
#   set.seed(54128)
  
  cat("Fitting model with orignal data...\n")
  m1 <- matrix_glm(x=x,Map=Map,formula=formula,family=family,response.name=response.name,verbose=verbose,...)
  
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
    
    m1 <- matrix_glm(x=Tab.bs,Map=Dat.bs,formula=formula,family=family,response.name=response.name,verbose=verbose,...)
    
    Res$boot$coefficients[,colnames(m1$coefficients),i] <- m1$coefficients
    Res$boot$SE[,colnames(m1$SE),i] <- m1$SE
  }
  
  Res$call <- match.call()
  Res$family <- family
  Res$pseudoreplicates <- SAMPLES
  Res$N <- N
  Res$coefficients <- apply(Res$boot$coefficients,c(1,2),mean,na.rm=TRUE)
  class(Res) <- "bootglm"
  
  return(Res)
}

bootstrap_glm.Dataset <- function(x,formula,family=poisson(link="log"),response.name="Count",
                                  N=100, verbose = FALSE...){
  m1 <- bootstrap_glm.default(x=x$Tab,Map=x$Map,formula=formula,family=family,response.name=response.name,N=N,...)
  m1$call <- match.call()
  return(m1)
}

summary.bootglm <- function(object,sortby="Variable",...){
  # Test data
#   object <- m1.bs
#   sortby <- "Variable"
  
  if(sortby != "Taxon" && sortby != "Variable"){
    stop("ERROR: You can only sort by Taxon or Variable",call.=TRUE)
  }
  
  TAB <- melt(object$coefficients,value.name="Estimate",varnames=c("Taxon","Variable"))
  quants <- apply(object$boot$coefficients,c(1,2),function(vector){x <- quantile(vector,probs=c(0.025,0.975),na.rm=TRUE)
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

print.summary.bootglm<-function(x,...){
  # Test
  #   x <- Res
  cat("Call:\n")
  print(x$call)
  
  row.names(x$coefficients) <- paste(as.character(x$coefficients$Taxon),as.character(x$coefficients$Variable),sep="__")
  x$coefficients$Taxon <- NULL
  x$coefficients$Variable <- NULL
  x$N <- format(x$N,digits=2)
  
  print(x$coefficients)
  
}


