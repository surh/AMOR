\name{read.pool}
\alias{read.pool}

\title{Read file describing which samples to pool}

\description{PROBABLY WILL BE DEPRECATED. Read file describing which samples to pool}

\usage{
read.pool(file)
}

\arguments{
  \item{file}{Path of text file describing which samples to pool. The file should be a tab delimited file, with each group as a line, and tab delimited entries with the ids of the samples to pool}
}
  
\details{}

\value{Returns a list with one entry per line of the file, and where each element is a vector containing the ids of the samples to pool}

\author{Sur from Dangl Lab}
