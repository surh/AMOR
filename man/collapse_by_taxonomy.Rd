\name{collapse_by_taxonomy}
\alias{collapse_by_taxonomy.default}
\alias{collapse_by_taxonomy.Dataset}

\title{Collapse taxa by taxonomy}

\description{Collapses the taxons in an abundance matrix according to their specified taxonomy}

\usage{
  \method{collapse_by_taxonomy}{default}(Tab,Tax,level=4,FUN=sum,sepchar=";")
  \method{collapse_by_taxonomy}{Dataset}(Dat,level=4,FUN=sum,sepchar=";")
}

\arguments{
  \item{Tab}{A \code{matrix} object encoding abundances of taxons (rows) on samples (columns).}
  \item{Tax}{A \code{data.frame} object encodeing the taxonomic information. Rows must be named and correspond to the names in Tab. The first to columns must be named "ID" and "Taxonomy". The column taxonomy must contain the taxonomic information of each taxon in the RDP format. See \link{create_dataset} for more info.}
  \item{Dat}{A \code{Dataset} object to be collapsed}.
  \item{level}{Taxonomic level to be used. Taxonomy must be in a format:
  <level1><sepchar><level2><sepchar><level3>...
  This option indicates which levels to use. Defaults to 4, which corresponds to phylum if the QIIME output from using the Greengenes database is used.}
  \item{FUN}{Function to apply when collapsing the data. Defaults to \code{sum}, which is ideal for count data. For proportioanl data \code{mean} or \code{median} migth be more appropriate. Any function that takes a vector of numbers and returns a single numeric value can be used.}
  \item{sepchar}{String that is used as separation of taxonomic levels. Defaults to the RDP standard ";"}
}

\details{Wrapper for \code{get_tax_level} and \code{collapse_matrix}.}

\value{The default method returns a \code{matrix} object.
The Dataset method returns a \code{Dataset} object when Dat includes a Map element (see \code{create_dataset}); when the Map element is missing it returns a \code{matrix} object.}

\references{RDP, Greengenes, QIIME}

\author{Sur from the Dangl Lab.}

\seealso{}

\examples{
library(AMOR)
data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)
Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)

# The following two are almost equivalent:
Tab.collapsed1 <- collapse_by_taxonomy(Dat$Tab,Dat$Tax)
Tab.collapsed2 <- collapse_by_taxonomy(Dat=Dat)
}
