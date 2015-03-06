\name{read.am}
\alias{read.am}

\title{Read abundance matrix}

\description{Reads multiple abundance matrix formats}

\usage{
read.am(file,format="am",taxonomy=FALSE)
}

\arguments{
  \item{file}{File path of the abundance matrix.}
  \item{format}{Format of the abundance table file. Either "am" or "qiime", see description for details.}
  \item{taxonomy}{Indicates whether a taxonomy assignment column is present in the abundance matrix. If FALSE, no taxonomy is expected, if a string or integer are passed, then the column with the corresponding name or number is expected to have taxonomy information.}
}
\details{
This function is a wrapper for \code{read.table()} reads tab delimited text files that contain abundance matrices. It supports two formats: 1) am, this format contains a first row of sample IDs and and a first colum of taxon IDs; and 2) qiime, This format is the original QIIME format, where the first line contains a version identification string, the second row sample IDs and the first column (after second row) taxon IDs. the first two lines start witht he '#' character.
}

\value{
Returns either a Dataset object, with Tax element NULL or extracted from the taxonomy column if option taxonomy is not FALSE.
}

\seealso{\code{\link{create_dataset}}}

\author{Sur from Dangl Lab}
