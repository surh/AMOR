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
