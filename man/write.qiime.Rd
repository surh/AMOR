\name{write.qiime}
\alias{write.qiime}
\alias{write.qiime.default}
\alias{write.qiime.Dataset}

\title{Write QIIME abundance matrix}

\description{Writes an abundance matrix in QIIME format.}

\usage{
\method{write.qiime}{default}(Tab, file)
\method{write.qiime}{Dataset}(Dat, file)
}
\arguments{
  \item{Tab}{Abundance matrix or data.frame.}
  \item{Dat}{A Dataset object.}
  \item{file}{Path for the file to create.}
}

\value{Returns the result of calling \code{\link{write.table()}}}

\author{Sur from Dangl Lab.}

