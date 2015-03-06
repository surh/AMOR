\name{clean}
\alias{clean}
\alias{remove_samples}
\alias{remove_taxons}

\title{Remove samples and taxons from a Dataset}

\description{Functions that will take a Dataset object, and remove samples and/or taxons with different criteria}

\usage{
clean(Dat,verbose=TRUE)
remove_samples(Dat,samples)
remove_taxons(Dat,taxons)
}

\arguments{
  \item{Dat}{A \code{Dataset} object.}
  \item{samples}{A character vector wiht the sample IDs to be removed.}
  \item{taxons}{A character vector with the taxon IDs to be removed.}
  \item{verbose}{logical indicating whether to print info about the results.}
}
\details{\code{clean} will remove any sample (column) or taxon (row) fromt your Dataset, that has no observations. Some functions like PCA and some linear models won't work if one of the rows or columns of the matrix are all zero, and so it is always a good idea to "clean" your dataset.

\code{remove_samples} and \code{remove_taxons} take respectively a list of samples or taxons to remove from the Dataset.
}

\value{Returns a \code{Dataset} object.}

\author{Sur from the Dangl Lab.}

\examples{
data(Rhizo)
data(Rhizo.map)
Dat <- create_dataset(Rhizo,Rhizo.map)
Dat <- remove_samples(Dat,colnames(Dat$Tab)[1:50])
Dat <- remove_taxons(Dat,row.names(Dat$Tab)[1:115])
sum(colSums(Dat$Tab) == 0)
sum(rowSums(Dat$Tab) == 0)
Dat$Tab
Dat.clean <- clean(Dat)
Dat.clean$Tab
}