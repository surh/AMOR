\name{pool_samples}
\alias{pool_samples.default}
\alias{pool_samples.Dataset}

\title{Pool samples}

\description{Pools abundances in samples according to grouping factor.}

\usage{
  \method{pool_samples}{default}(Tab,groups,FUN=sum)
  \method{pool_samples}{Dataset}(Dat,groups,FUN=sum)
}

\arguments{
  \item{Tab}{A \code{matrix} object encoding abundances of taxons (rows) on samples (columns).}
  \item{Dat}{A \code{Dataset} object to be collapsed}.
  \item{groups}{For the default method this is a vector or factor specifying to which group each sample belongs. Vectors will be converted to factorrs with the \code{factor} function.
  For the Dataset method, this is a single character or numeric value indicating the column from Dat$Map to be used as a grouping factor.}
  \item{FUN}{Function to apply when collapsing the data. Defaults to \code{sum}, which is ideal for count data. For proportioanl data \code{mean} or \code{median} migth be more appropriate. Any function that takes a vector of numbers and returns a single numeric value can be used.}
}

\details{Wrapper for \code{collapse_matrix}.}

\value{The default method returns a \code{matrix} object.
The Dataset method returns a \code{Dataset} object when Dat includes a Tax element (see \code{create_dataset}); when the Map element is missing it returns a \code{matrix} object.}

\references{}

\author{Sur from the Dangl Lab.}

\seealso{}

\examples{
library(AMOR)
data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)
Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)

# The following are almost equivalent.
Tab.collapsed1 <- pool_samples(Tab = Dat$Tab,groups = Dat$Map$fraction)
Tab.collapsed2 <- pool_samples(Dat = Dat,groups = "fraction")
}
