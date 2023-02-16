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
#' @details Implments abundance filters defined in Benson et al.(2010) & Lundberg et al. (2012).
#' It is recommended to threshold any Dataset before any quantitative analysis.
#' 
#' @return When table = FALSE, it returns a logical index vector that shows which taxa pass the
#' thresholds.
#' 
#' When table = TRUE, it returns a matrix or Dataset object (depending on the method used),
#' where samples that do not pass the thresholds have been removed
#' 
#' @author Sur from Dangl Lab.
#' 
#' @references 
#' 1. Benson AK, Kelly S a, Legge R, Ma F, Low SJ, Kim J, et al. Individuality in gut microbiota
#' composition is a complex polygenic trait shaped by multiple environmental and host genetic
#' factors. Proc Natl Acad Sci U S A. 2010 Nov 2;107(44):18933–8. Available from:
#' \url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2973891&tool=pmcentrez&rendertype=abstract}
#' 
#' 2. Lundberg DS, Lebeis SL, Herrera Paredes S, Yourstone S, Gehring J, Malfatti S, et al.
#' Defining the core Arabidopsis thaliana root microbiome. Nature. Nature Publishing Group;
#' 2012 Aug 1;488(7409):86–90. Available from: \url{http://www.nature.com/doifinder/10.1038/nature11237}
#' 
#' @export
#' @name measurable_taxa
measurable_taxa <- function(...) UseMethod("measurable_taxa")

#' @rdname measurable_taxa
#' @name measurable_taxa
#' @method measurable_taxa default
#' @export
measurable_taxa.default <- function(Tab, min_reads_otu = 25, min_samples_otu = 5,
                                    method = "absolute", table = TRUE){
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
  
  if(sum(index) == 0){
    cat("No taxa passes the specified threshold.\n")
    return(NULL)
  }
  
  if(table){
    Tab <- Tab[ index, , drop = FALSE ]
    return(Tab)
  }else{
    return(index)
  }
}

#' @rdname measurable_taxa
#' @name measurable_taxa
#' @method measurable_taxa Dataset
#' @export
measurable_taxa.Dataset <- function(Dat, min_reads_otu = 25, min_samples_otu = 5, 
                                    method="absolute", table = TRUE,clean = TRUE){
  
  res <- measurable_taxa(Tab = Dat$Tab, min_reads_otu = min_reads_otu,
                         min_samples_otu = min_samples_otu, method = method,
                         table = FALSE)
  
  if(is.null(res)){
    return(res)
  }
  
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

#' @export
findGoodOTUs <- function(...){
  warning("This function will be deprecated soon. Use measurable_taxa instead\n",call. = TRUE)
  measurable_taxa(...)
}

