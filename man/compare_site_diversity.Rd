\name{compare_site_diversity}
\alias{compare_site_diversity}
\alias{compare_site_diversity.default}
\alias{compare_site_diversity.Dataset}

\title{Compare diversity accross sites and groups of samples.}

\description{PErforms permuation to estimate diversity as a function of adding sites. And does this for samples grouped according to a factor.}

\usage{
  \method{compare_site_diversity}{default}(Tab, Map, factor,
    divfun = total_richness, nperm = 20)
  \method{compare_site_diversity}{Dataset}(Dat, factor,
    divfun = total_richness, nperm = 20)
}

\arguments{
  \item{Tab}{A numeric matrix of samples as columns and taxa as rows.}
  \item{Map}{A data.frame containing the variables to be modelled as columns and samples as rows. The rows should be named with sample IDs and must correspond to the column names from x if an abundance matrix was passed}
  \item{Dat}{A Dataset object}
  \item{factor}{String representing the name of the variable to be used for grouping samples. Must correspond to a header name in the \code{Mat} portion of the Dataset object.}
  \item{divfun}{Function that returns a diversity estimate given a matrix of samples. See \link{total_richness} and examples to see how to define your function.}
  \item{nperm}{Number of permutations to perform.}

}

\value{A data.frame of class site.diversity which contains the following variables:
  \item{mean}{Mean diversity value of all permutations}
  \item{sd}{Standard deviation of the diversity estimates estimated from the permutations.}
  \item{nsites}{Number of sites (ie. samples).}
  \item{group}{Variable used for selecting the samples.}
}

\author{Sur from Dangl Lab}

\seealso{\link{site_diversity}, \link{total_richness}}

\examples{
data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)
Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)

sitediv.accession <- compare_site_diversity(Dat = Dat,factor = "accession", divfun = total_richness, 20)
plotgg(sitediv.accession)

# The following code requires the vegan package
#divfun <- function(x){
#  if(!is.null(ncol(x)))
#    x <- rowSums(x)
#  s <- vegan::diversity(x)
#  return(s)
#}

#sitediv.accession <- compare_site_diversity(Dat = Dat,factor = "accession", divfun = divfun, 20)
#plotgg(sitediv.accession, alpha = 0.3) + 
#  scale_color_brewer(palette = "Set3") +
#  scale_fill_brewer(palette = "Set3")
}
