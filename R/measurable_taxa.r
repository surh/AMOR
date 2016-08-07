#' Identify measurable taxa
#' 
#' Identfies taxa above certain abundance and prevalence
#' thresholds, that can be quantitatively analyzed
#' 
#' @aliases findGoodOTUs
#' 
#' @param Tab Count table matrix. Samples as columns and OTUs as rows.
#' @param Dat A \code{Dataset} object.
#' @param method Either "absolute" (min absolute number of reads in min
#' number of sample) or "AbsProp" (absolute number if reads and proportion
#' of samples).
#' @param table logical indicating whether an abundance table should be
#' returned
#' @param clean logical, indicating whether the function \code{\link{clean}}
#' should be applied to the result.
#' 
#' @return When table = FALSE, it returns a logical index vector that shows which taxa pass the
#' thresholds.
#' 
#' When table = TRUE, it returns a matrix or Dataset object (depending on the method used),
#' where samples that do not pass the thresholds have been removed
#' 
#' @author Sur from Dangl Lab.
#' 
#' @export
measurable_taxa <- function(...) UseMethod("measurable_taxa")

#' @rdname measurable_taxa
#' @method measurable_taxa default
measurable_taxa.default <- function(Tab,min_reads_otu,min_samples_otu,method="absolute",
                              table = TRUE){
  if(method == "absolute"){
    #Count samples where OTU is above threshold
    samples_per_otu <- rowSums(Tab >= min_reads_otu)
    index <- samples_per_otu >= min_samples_otu
    #return(index)
  }else if(method == "AbsProp"){
    if(min_samples_otu < 0 || min_samples_otu > 1){
      stop("findGoodOTUs: Invalid proportion of samples. Must be in the range [0,1]",call.=TRUE)
    }
    reads_per_otu <- rowSums(Tab)
    #Count samples in which OTU appears
    samples_per_otu <- rowSums(Tab > 0)
    #Convert proportion of samples to absolute number
    min_samples_otu <- min_samples_otu*dim(Tab)[2]
    print(min_samples_otu)
    index1 <- samples_per_otu >= min_samples_otu
    index2 <- reads_per_otu >= min_reads_otu
    index <- index1 & index2
    #return(index)
  }else{
    stop("findGoodOTUs: Invalid method",call.=TRUE)
  }  
  
  if(table){
    Tab <- Tab[ index, , drop = FALSE ]
    return(Tab)
  }else{
    return(index)
  }
}

#' @rdname measurable_taxa
#' @method measurable_taxa Dataset
measurable_taxa.Dataset <- function(Dat,min_reads_otu,min_samples_otu,method="absolute",
                        table = TRUE,clean = TRUE){
  
  
  # Dat = Dat
  # min_reads_otu = 25
  # min_samples_otu = 8
  # table <- TRUE
  # clean <- TRUE
  # method <- "absolute"
  
  res <- measurable_taxa(Tab = Dat$Tab, min_reads_otu = min_reads_otu,
                         min_samples_otu = min_samples_otu, method = method,
                         table = FALSE)
  
  if(table){
    res <- create_dataset(Tab = Dat$Tab[ res, , drop = FALSE ],
                          Map = Dat$Map,
                          Tax = Dat$Tax[ res, , drop = FALSE] )
    if(clean){
      res <- clean(res)
    }
  }
  
  return(res)
}

findGoodOTUs <- function(...){
  warning("This function will be deprecated soon. Use measurable_taxa instead\n",call. = TRUE)
  measurable_taxa(...)
}

