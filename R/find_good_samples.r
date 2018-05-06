#' Find good samples
#' @aliases findGoodSamples
#' 
#' Find samples with minumum number of reads
#'
#' @param x Abundance matrix with samples as
#' columns
#' @param min_counts Minimum count per sample to keep
#'
#' @return A boolean vector specifying which samples to keep
#' 
#' @export
#' @author Sur Herrera Paredes
#'
#' @examples
#' mat <- matrix(rep(c(2,0,1),each = 10), ncol = 3)
#' mat[ , find_good_samples(mat, 1)]
find_good_samples <- function(x,min_counts){
  # Finds OTUs with reads above threshold
  index <- colSums(x) > min_counts
  return(index)
}

findGoodSamples <- find_good_samples