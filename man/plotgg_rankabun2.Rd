\name{plotgg_rankabun2}
\alias{plotgg_rankabun2}
\alias{plotgg_rankabun2.default}
\alias{plotgg_rankabun2.Dataset}


\title{Make rank abundance plot with Poisson confidence intervals.}

\description{Make a rank abundance plot with Poisson confidence intervals from a Dataset object.

IMPORTANT: There is no guarantee that Poisson confidence intervals are representative of the data. This should be viewed only as descriptive statistic and not for inference.}

\usage{
  \method{plotgg_rankabun2}{default}(Tab, Map, groupby, sortby,
                                    alpha = 0.2, theme = theme_blackbox, 
                                    variable.name = "Taxon", value.name = "Abundance", 
                                    sample.id.name = "SAMPLEID")                                    
  \method{plotgg_rankabun2}{Dataset}(Dat, groupby, sortby,
                                    alpha = 0.2, theme = theme_blackbox, 
                                    variable.name = "Taxon", value.name = "Abundance", 
                                    sample.id.name = "SAMPLEID")

}

\arguments{
  \item{Tab}{A matrix object with samples as columns and taxa as rows.}
  \item{Map}{A mapping file for \code{Tab}. Each row must be a sample with row names matching column names in \code{Tab} and in the same order as in \code{Tab}.}
  \item{Dat}{A Dataset object.}
  \item{groupby}{Variable name to be used for grouping samples before plotting th rank abundance.}
  \item{sortby}{Variable value to be used as reference for sorting the taxa in the rank abundance plot.}
  \item{alpha}{Transparency parameter for ggplot2}
  \item{theme}{ggplot2 theme to be used for plotting.}
  \item{variable.name}{x-axis label in the plot.}
  \item{value.name}{y-axis label in the plot.}
  \item{FUN}{Function to use to aggregate samples according to \code{groupby}}
  \item{sample.id.name}{A name to store sample IDs in \code{Map}. Used for internal handling only.}
}
\details{}

\value{A ggplot2 plot.}

\author{Sur from Dangl Lab}

\seealso{\link{plotgg_rankabun}}

\examples{
data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)
Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)

plotgg_rankabun2(Tab = Dat$Tab,Map = Dat$Map,groupby = "fraction",sortby = "E")
plotgg_rankabun2(Dat,groupby = "fraction",sortby = "E")
}
