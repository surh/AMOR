\name{rarefaction}
\alias{rarefaction}

\title{Rarefaction}

\description{Performs rarefaction on count table.}

\usage{
  \method{rarefaction}{default}(x, sample)
  \method{rarefaction}{Dataset}(x,sample)
}

\arguments{
  \item{x}{Either a matrix object or a \code{Dataset} object. If a matrix is passed, samples must be the columns and rows the species.}
  \item{sample}{Sample size, either a single integer value or a vector of the same length the number of columns (samples) of x.}
}

\details{Modified from the rrarefy function of vegan-2.0-4.}

\value{The defualt method returns a matrix with the rarefied counts. The Dataset method returns a \code{Dataset } object where the Tab element is the rarefied counts.}

\author{Sur from Dangl Lab.}

\examples{
data(Rhizo)
data(Rhizo.map)
Dat <- create_dataset(Rhizo,Rhizo.map)
set.seed(712) # Always specify and save the seed before rarefaction
rarefaction(x=Dat$Tab,sample=100)
set.seed(712)
rarefaction(x=Dat,sample=100)$Tab
}
