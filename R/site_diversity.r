#' Diversity across sites.
#' 
#' Calculates a diversity measure as a function of adding samples,
#' and uses permutation to obtain a confidence interval.
#'
#' @param Dat A Dataset object
#' @param factor String representing the name of the variable to be
#' used for grouping samples. Must correspond to a header name in
#' the Mat portion of the Dataset object.
#' @param group Value of the variable defined in factor to extract.
#' See examples.
#' @param divfun Function that returns a diversity estimate given
#' a matrix of samples. See total_richness and examples to see how
#' to define your function.
#' @param nperm Number of permutations to perform.
#'
#' @return A data.frame of class \code{site.diversity} which contains the following variables:
#' \describe{
#'   \item{mean}{Mean diversity value of all permutations}
#'   \item{sd}{Standard deviation of the diversity estimates estimated from the permutations.}
#'   \item{nsites}{Number of sites (ie. samples).}
#'   \item{group}{Variable used for selecting the samples.}
#' }
#' @export
#' @author Sur Herrera Paredes
#' @seealso \link{compare_site_diversity}, \link{total_richness}
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' data(Rhizo.tax)
#' Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)
#' 
#' soil.sitediv <- site_diversity(Dat = Dat,
#'                                factor = "accession",
#'                                group = "Soil")
#' col.sitediv <- site_diversity(Dat = Dat,
#'                               factor = "accession",
#'                               group = "Col")
#' plotgg(soil.sitediv)
#' plotgg(col.sitediv)
#' 
#' # The following function requires the vegan package installed
#' # and can be used with the site_diversity function.
#' #divfun <- function(x){
#' #  if(!is.null(ncol(x)))
#' #    x <- rowSums(x)
#' #  s <- vegan::diversity(x)
#' #  return(s)
#' #}
site_diversity <- function(Dat,factor,group,divfun = total_richness,nperm=20){
  
  divfun <- match.fun(divfun)
  Dat.temp <- remove_samples(Dat, row.names(Dat$Map)[ Dat$Map[,factor] != group ])
  #Dat.temp$Tab <- Dat.temp$Tab > 0
  Perms <- matrix(ncol = nrow(Dat.temp$Map), nrow = nperm)
  for (i in 1:nperm){
    #i <- 1
    permutation <- sample(row.names(Dat.temp$Map))
    for (j in 1:length(permutation)){
      #j <- 2
      mat <- matrix(Dat.temp$Tab[,permutation[1:j]],ncol = j)
      S <- divfun(mat)
      Perms[i,j] <- S
    }
    
  }
  #Perms
  Res <- apply(Perms,2,function(vec) data.frame(mean = mean(vec), sd = sd(vec)))
  Res <- do.call(rbind,Res)
  Res$nsites <- 1:j
  Res$group <- group
  class(Res) <- c("site.diversity","data.frame")
  return(Res)
}

#' Compare diversity accross sites and groups of samples.
#' 
#' Performs permuation to estimate diversity as a function
#' of adding sites. And does this for samples grouped
#' according to a factor.
#'
#' @param Tab A numeric matrix of samples as columns and taxa as rows.
#' @param Map A data.frame containing the variables to be modelled as
#' columns and samples as rows. The rows should be named with sample
#' IDs and must correspond to the column names from x if an abundance
#' matrix was passed
#' @param Dat A Dataset object.
#' @param factor String representing the name of the variable to be
#' used for grouping samples. Must correspond to a header name in the
#' Mat portion of the Dataset object.
#' @param divfun Function that returns a diversity estimate given a
#' matrix of samples. See total_richness and examples to see how to
#' define your function.
#' @param nperm Number of permutations to perform.
#'
#' @return A data.frame of class \code{site.diversity} which contains the following variables:
#' \describe{
#'   \item{mean}{Mean diversity value of all permutations}
#'   \item{sd}{Standard deviation of the diversity estimates estimated from the permutations.}
#'   \item{nsites}{Number of sites (ie. samples).}
#'   \item{group}{Variable used for selecting the samples.}
#' }
#' @export
#' @author Sur Herrera Paredes
#' @seealso \link{site_diversity}, \link{total_richness}
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' data(Rhizo.tax)
#' Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)
#' 
#' sitediv.accession <- compare_site_diversity(Dat = Dat,factor = "accession", divfun = total_richness, 20)
#' plotgg(sitediv.accession)
#' 
#' # The following code requires the vegan package
#' # divfun <- function(x){
#' #  if(!is.null(ncol(x)))
#' #    x <- rowSums(x)
#' #  s <- vegan::diversity(x)
#' #  return(s)
#' # }
#' 
#' # sitediv.accession <- compare_site_diversity(Dat = Dat,
#' #                                             factor = "accession",
#' #                                             divfun = divfun, 20)
#' # plotgg(sitediv.accession, alpha = 0.3) +
#' #   scale_color_brewer(palette = "Set3") +
#' #   scale_fill_brewer(palette = "Set3")
compare_site_diversity <- function(...) UseMethod("compare_site_diversity")

#' @rdname compare_site_diversity
#' @method compare_site_diversity default
compare_site_diversity.default <- function(Tab, Map, factor,
                                           divfun = total_richness,
                                           nperm = 20){
  
  groups <- levels(Map[,factor])
  Dat <- create_dataset(Tab = Tab, Map = Map)
  Res <- NULL
  for(group in groups){
    res <- site_diversity(Dat = Dat, factor = factor,
                          group = group, divfun = divfun,
                          nperm = nperm)
    Res <- rbind(Res,res)
  }
  
  class(Res) <- c("site.diversity","data.frame")
  return(Res)
}

#' @rdname compare_site_diversity
#' @method compare_site_diversity Dataset
compare_site_diversity.Dataset <- function(Dat, factor,
                                           divfun = total_richness,
                                           nperm = 20){
  
  Res <- compare_site_diversity(Tab = Dat$Tab,
                                Map = Dat$Map,
                                factor = factor,
                                divfun = divfun,
                                nperm = nperm)
  Res$group <- factor(Res$group, levels = levels(Dat$Map[,factor]))
  return(Res)
}