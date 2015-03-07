\name{PCA}
\alias{PCA}
\alias{PCA.default}
\alias{PCA.Dataset}


\title{Principal Component Analysis}

\description{Function that performs principal component analysis on an abundance matrix}

\usage{
  \method{PCA}{default}(x, cor = FALSE, dim = min(nrow(x), ncol(x)))
  \method{PCA}{Dataset}(Dat, cor = FALSE, dim = min(nrow(Dat$Tab), ncol(Dat$Tab)))
}

\arguments{
  \item{x}{Numeric matrix where samples are columns and rows are species.}
  \item{cor}{logical value indicating whether the correlation matrix should be used instead of the covariance matrix.}
  \item{dim}{Number of dimensions to return.}
  \item{Dat}{A Dataset object, see \code{\link{create_dataset}}.}
}

\details{This function is the same as function \code{\link{pca}} from the \code{\link{labdsv}} package, but includes a methdod for Dataset objects.}

\value{A \code{PCA} object. Includes the same attributes as a \code{\link{pca}} object from the \code{\link{labdsv}} package. When the Dataset method is used, it includes two additional slots:
  \item{Map}{The Mapping file for the samples.}
  \item{Tax}{The Taxonomic information of the taxa.}
}

\author{Sur from Dangl Lab.}

\seealso{\link{create_dataset}, \link{pca},\link{PCO},\link{pco},\link{plotgg.pca}}

\examples{
data(Rhizo)
data(Rhizo.map)
Dat <- create_dataset(Rhizo,Rhizo.map)
Dat.pca <- PCA(Dat)
plotgg(Dat.pca,col="accession",shape="fraction",point_size=4,biplot=TRUE)
}

