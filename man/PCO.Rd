\name{PCO}
\alias{PCO}
\alias{PCO.default}
\alias{PCO.Dataset}


\title{Principal Coordinates Analysis}

\description{Function that performs Principal Coordinates Analysis on abundance matrix.}

\usage{
  \method{PCO}{default}(x, dim = 3)
  \method{PCO}{Dataset}(Dat, dim = 3, distfun = dist)
}

\arguments{
  \item{x}{Distance matrix, must be a \code{dist} object.}
  \item{dim}{Number of dimensions to return.}
  \item{Dat}{\code{Dataset} object. See \code{\link{create_dataset}} for more information.}
  \item{distfun}{Function that calculates distance matrix for abundance table in \code{Dataset} object.}
}

\details{This is the same as function \code{\link{pco}} from the labdsv, but it includes a \code{Dataset} method.}

\value{A PCO object as defined by the \code{\link{pco}} function.}

\author{Sur from Dangl Lab.}

\seealso{
  \code{\link{create_dataset}}, \code{\link{pca}},\code{\link{PCA}},\code{\link{pco}},\code{\link{dsvdis}}
}

\examples{
data(Rhizo)
data(Rhizo.map)
Dat <- create_dataset(Rhizo,Rhizo.map)

# distfun <- function(x) vegan::vegdist(x,method="bray") #requires vegan package
distfun <- dist

Dat.pco <- PCO(Dat,dim=2,distfun=distfun)
plotgg(Dat.pco)
plotgg(Dat.pco,shape="fraction",point_size=3)
plotgg(Dat.pco,shape="fraction",col="accession",point_size=4)

}


