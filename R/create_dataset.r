create_dataset <- function(Tab=NULL,Map=NULL,Tax=NULL){
  # Test data
#   Tab <- Rhizo
#   Map <- Rhizo.map
#   Tax <- Rhizo.tax
  
  Dataset <- list(Tab = NULL,Map = NULL,Tax = NULL)
  
  # Make sure Tab is a numeric matris
  Tab <- as.matrix(Tab)
  if(!is.numeric(Tab)){
    stop("ERROR: Tab is not, or cannot be coerced into a numeric matrix",call.=TRUE)
  }
  Dataset$Tab <- Tab
  
  # Make sure Map and Tab are consistent, and add Map if exists
  if(!is.null(Map)){
    Map <- as.data.frame(Map)
    if(any(colnames(Tab) != row.names(Map)) || ncol(Tab) != nrow(Map)){
      stop("ERROR: Samples in Tab and Map do not match",call.=TRUE)
    }
    Dataset$Map <- Map
  }
    
  # If Tax is passed add it
  if(!is.null(Tax)){
    # Make sure Tab and Tax are consistent
    Tax <- as.data.frame(Tax)
    if(any(row.names(Tab) != row.names(Tax)) || nrow(Tab) != nrow(Tax)){
      stop("ERROR: Taxons in Tab and Tax do not match",call.=TRUE)
    }
    Dataset$Tax <- Tax
  }
  
  class(Dataset) <- "Dataset"
  return(Dataset)
}

#### UTILS ####
print.Dataset <- function(x){
  if(!any("Dataset" %in% class(x))){
    stop("ERROR: Must pass a Dataset object.",call. = TRUE)
  }
  cat("A Dataset object:\n")
  
  nsamples <- length(samples(x))
  ntaxa <- length(taxa(x))
  nvars <- length(variables(x))
  
  cat("  There are ", ntaxa, " taxa in ",nsamples," samples\n",sep = "")
  if(nvars > 0){
    cat("  There are ", nvars, " metadata variables\n",sep = "")
  }
  
  cat("\n  Samples: ",paste(samples(x)[1:min(5,nsamples)],collapse = ","),"...\n",sep = "")
  cat("  Taxa: ",paste(taxa(x)[1:min(5,ntaxa)],collapse = ","),"...\n", sep = "")
  if(nvars > 0){
    cat("  Variables: ",paste(variables(x)[1:min(5,nvars)],collapse = ","),"...\n", sep = "")
  }
  
  # Most relevant otu
  mos.abun <- which.max(rowSums(x$Tab))
  mos.prev <- which.max(rowSums(x$Tab > 0))
  
  cat("\n  The most abundant taxon is (are): ", paste(taxa(x)[ mos.abun ], collapse = ","),"\n",sep = "")
  cat("  The highet abundance is: ", sum(x$Tab[mos.abun[1],]),"\n",sep = "")
  
  cat("\n  The most prevalent taxon is (are): ", paste(taxa(x)[ mos.prev ], collapse = ","),"\n",sep = "")
  cat("  The highet prevalence is: ", sum(x$Tab[mos.prev[1],] > 0),sep = "")
}

#' Get sample names
#' 
#' Get names of samples in Dataset
#' 
#' @param Dat A dataset object
#' 
#' @return A vector of sample names
#' @seealso \code{\link{create_dataset}}, \code{\link{taxa}},
#' \code{\link{variables}}
samples <- function(Dat){
  if(class(Dat) != "Dataset"){
    stop("ERROR: You must pass a Dataset object",call. = TRUE)
  }
  samples <- colnames(Dat$Tab)
  return(samples)
}

#' Get metadata variables
#' 
#' Get names of variables in metadata
#' 
#' @param Dat A dataset object
#' 
#' @return A vector of variable names
#' 
#' @seealso \code{\link{create_dataset}}, \code{\link{taxa}},
#' \code{\link{samples}}
variables <- function(Dat){
  if(class(Dat) != "Dataset"){
    stop("ERROR: You must pass a Dataset object",call. = TRUE)
  }
  vars <- colnames(Dat$Map)
  return(vars)
}


#' Get taxa names
#' 
#' Get names of taxa in Dataset
#' 
#' @param Dat A dataset object
#' 
#' @return A vector of taxa names
#' @seealso \code{\link{create_dataset}}, \code{\link{variables}},
#' \code{\link{samples}}
taxa <- function(Dat){
  if(class(Dat) != "Dataset"){
    stop("ERROR: You must pass a Dataset object",call. = TRUE)
  }
  taxa <- row.names(Dat$Tab)
  return(taxa)
}

