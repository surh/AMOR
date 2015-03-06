\name{plotgg_taxon}
\alias{plotgg_taxon}
\alias{plotgg_taxon.default}
\alias{plotgg_taxon.Dataset}

\title{Plot single taxon.}

\description{Make a box and whiskers plot of a single taxon, while grouping by a variable.}

\usage{
  \method{plotgg_taxon}{default}(Tab, Map, taxon ,x, col = NULL,
                        var.name = "Abundance",theme = theme_blackbox)
  \method{plotgg_taxon}{Dataset}(Dat, taxon ,x, col = NULL,
                        var.name = "Abundance",theme = theme_blackbox)
}

\arguments{
  \item{Tab}{Matrix with samples as columns and taxa as rows.}
  \item{Map}{Mapping file for \code{Tab}. Each row must correspond to a sample in \code{Tab}, they must be in the same order, and be named with the same names as the columns in \code{Tab}}
  \item{Dat}{A Dataset object.}
  \item{taxon}{String indicating which taxon to plot.}
  \item{x}{String indicating which variable to use as the x-axis in the plot. Should correspond to a header name in \code{Map}}
  \item{col}{String indicating which variable to use to color the plot. Can be the same or different as \code{x}. It should correspond to a header name in \code{Map}.}
  \item{var.name}{y-axis label on the plot.}
  \item{theme}{A ggplot2 theme to be used woth the plot.}
}

\details{}

\value{A ggplot2 plot.}

\author{Sur from Dangl Lab}

\examples{
data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)
Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)

# Passing independently an abundance table and mapping file
plotgg_taxon(Tab = Dat$Tab, Map = Dat$Map, taxon = "OTU_14834", x = "fraction")
plotgg_taxon(Tab = Dat$Tab, Map = Dat$Map, taxon = "OTU_14834", x = "fraction", col = "accession")

# Passing a dataset object
plotgg_taxon(Dat, taxon = "OTU_14834", x = "fraction")
plotgg_taxon(Dat, taxon = "OTU_14834", x = "fraction", col = "accession")
}
