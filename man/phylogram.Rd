\name{phylogram}
\alias{phylogram}
\alias{phylogram.default}
\alias{phylogram.Dataset}

\title{Phylogram}

\description{Make a relative abundance barplot of taxon abundances.}

\usage{
  \method{phylogram}{default}(Tab,Map=NULL,facet = NULL,colname="Sample",
  variable.name="Taxon",value.name="Abundance",
  scales="free_x",nrow.legend=20,ntaxa = NULL,
  other_name = "other")
  \method{phylogram}{Dataset}(Dat,facet = NULL,colname="Sample",
                              variable.name="Taxon",value.name="Abundance",
                              scales="free_x",nrow.legend=20,ntaxa = NULL,
                              other_name = "other")
}

\arguments{
  \item{Tab}{A matrix object representing samples as columns and taxa as rows.}
  \item{Dat}{A \code{Dataset} object}
  \item{Map}{A mapping file matching the \code{Tab} parameter. One row per sample, and the rows must be named to match the column names in \code{Tab}, and they must be in the same order.}
  \item{facet}{facet formula for \code{\link{facet_grid}}.}
  \item{colname}{Label to be used in the x-axis of the phylogram.}
  \item{variable.name}{Label to be used in for the colors in the phylogram.}
  \item{value.name}{Label for the y-axis of the phylogram.}
  \item{scales}{scales option for \code{\link{facet_grid}}.}
  \item{nrow.legend}{number of rows to use in the legend.}
  \item{ntaxa}{Number of taxa to plot, the rest will be collapsed into one category. Samples will be sorted by abundance and the top \code{ntaxa} will be plotted.}
  \item{other_name}{Name to give to the collapsed taxa when \code{ntaxa} is used.}
}

\value{ A ggplot2 object of the plot.}

\author{Sur from Dangl Lab}

\examples{
data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)
Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)

Dat.phyl <- collapse_by_taxonomy(Dat,4)
phylogram(Tab = Dat.phyl$Tab,Map = Dat.phyl$Map,facet = ~ fraction)
phylogram(Dat.phyl,facet = ~ fraction)
phylogram(Dat.phyl,facet = ~ fraction, ntaxa = 5)
}

