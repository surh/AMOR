% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.am.r
\name{read.am}
\alias{read.am}
\title{Read abundance matrix}
\usage{
read.am(file, format = "am", taxonomy = FALSE, simplify = FALSE)
}
\arguments{
\item{file}{File path of the abundance matrix.}

\item{format}{Format of the abundance table file.
Either "am" or "qiime", see description for details.}

\item{taxonomy}{Indicates whether a taxonomy assignment
column is present in the abundance matrix. If FALSE,
no taxonomy is expected, if a string or integer are passed,
then the column with the corresponding name or number is
expected to have taxonomy information.}

\item{simplify}{Logical, indicating whether the abundance
table should be coerced to a numerical matrix}
}
\value{
If simplify = FALSE, it returns a Dataset object, where the value
of the Tax element is determined by the taxonomy argiment.

If simplify = TRUE, it returns a numeric abundance matrix. Any taxonomic
information is lost in this case.
}
\description{
Reads an abundance matrix into a Dataset object
}
\details{
This function is a wrapper for read.table(), and
reads tab-delimited text files that contain abundance matrices.
It supports two formats: 1) am, this format contains a first row
of sample IDs and and a first colum of taxon IDs; and 2) qiime,
this format is the original QIIME format, where the first line
contains a version identification string, the second row sample
IDs and the first column (after second row) taxon IDs. The first
two lines start witht he '#' character.
}
\seealso{
\code{\link{create_dataset}}
}
\author{
Sur Herrera Paredes
}
