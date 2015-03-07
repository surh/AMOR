\name{total_richness}
\alias{total_richness}

\title{Calculate richness of an abundance table.}

\description{Given a numeric matri representign an abundance table, it calculates the total richness (ie species number) in it.}

\usage{
  total_richness(Tab)
}

\arguments{
  \item{Tab}{A numeric matrix with columns as samples and rows as taxa.}
}

\value{A number representing the total richness.}

\author{Sur from Dangl Lab}

\seealso{\link{site_diversity}, \link{compare_site_diversity}}
