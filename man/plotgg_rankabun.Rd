\name{plotgg_rankabun}
\alias{plotgg_rankabun}
\alias{plotgg_rankabun.default}
\alias{plotgg_rankabun.Dataset}

\title{Make rank abundance plot}

\description{Make a rank abundance plot from a Dataset object.}

\usage{
  \method{plotgg_rankabun}{default}(Tab, Map = NULL, groupby = NULL, sortby = NULL,
                                    theme = theme_blackbox, variable.name = "Taxon",
                                    value.name = "Abundance", FUN = mean)
  \method{plotgg_rankabun}{Dataset}(Dat, groupby = NULL, sortby = NULL,
                                    theme = theme_blackbox, variable.name = "Taxon",
                                    value.name = "Abundance", FUN = mean)

}

\arguments{
  \item{Tab}{A matrix object with samples as columns and taxa as rows.}
  \item{Map}{A mapping file for \code{Tab}. Each row must be a sample with row names matching column names in \code{Tab} and in the same order as in \code{Tab}.}
  \item{Dat}{A Dataset object.}
  \item{groupby}{Variable name to be used for grouping samples before plotting th rank abundance.}
  \item{sortby}{Variable value to be used as reference for sorting the taxa in the rank abundance plot.}
  \item{theme}{ggplot2 theme to be used for plotting.}
  \item{variable.name}{x-axis label in the plot.}
  \item{value.name}{y-axis label in the plot.}
  \item{FUN}{Function to use to aggregate samples according to \code{groupby}}
}
\details{}

\value{A ggplot2 plot.}

\author{Sur from Dangl Lab}

\seealso{\link{plotgg_rankabun2}}

\examples{
data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)
Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)

plotgg_rankabun(Tab = Dat$Tab)
plotgg_rankabun(Tab = Dat$Tab, Map = Dat$Map,groupby = "fraction",sortby = "Soil")

plotgg_rankabun(Dat)
plotgg_rankabun(Dat, groupby = "fraction", sortby = "E")
}
