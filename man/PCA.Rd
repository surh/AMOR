\name{PCA}
\alias{PCA}

\title{Principal Component Analysis}

\description{Function that performs principal component analysis on an abundance matrix}

\usage{
  \method{PCA}{default}(mat, cor = FALSE, dim = min(nrow(mat), ncol(mat)))
  \method{PCA}{Dataset}(Dat,...)
}

\arguments{
  \item{mat}{Matrix where samples are columns and rows are species.}
  \item{cor}{logical value indicating whether the correlation matrix should be used instead of the covariance matrix.}
  \item{dim}{Number of dimensions to return.}
  \item{Dat}{A Dataset object, see \code{\link{create_dataset}}.}
  \item{dim}{Other arguments to pass to the default method.}
}

\details{This function is the same as function \code{\link{pca}} from the \code{\link{labdsv}} package, but includes a methdod for Dataset objects.}

\value{A \code{PCA} object. Includes the same attributes as a \code{\link{pca}} object from the \code{\link{labdsv}} package. When the Dataset method is used, it includes two additional slots:
  \item{Map}{The Mapping file for the samples.}
  \item{Tax}{The Taxonomic information of the taxa.}
}

\author{Sur from Dangl Lab.}

\seealso{
  \code{\link{create_dataset}}, \code{\link{pca}},\code{\link{PCO}},\code{\link{pco},\code{\link{plotgg.pca}}}
}

\examples{
data(Rhizo)
data(Rhizo.map)
Dat <- create_dataset(Rhizo,Rhizo.map)
Dat.pca <- PCA(Dat)
plotgg(Dat.pca,col="accession",shape="fraction",point_size=4,biplot=TRUE)
}

