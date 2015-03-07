\name{heatgg}
\alias{heatgg}
\alias{heatgg.default}
\alias{heatgg.Dataset}

\title{Heatmap with ggplot2%%  ~~function to do ... ~~
}
\description{Basic heatmap function with ggplot2. I allows either to sort the samples accordin to a user defined variable, or do a hierarchichal clustering and uses ggdendro and grid to plot it}

\usage{
\method{heatgg}{default}(Tab, Map, order.samples.by = NULL, facet = NULL,
      sample.id.var = "SAMPLEID", col.name = "Taxon",
      value.name = "Abundance", trans = "log", guide = "colourbar",
      gradientn.colours = c("white","#67001F"), facet.scales = "free",
      cluster = FALSE, distfun = dist )
\method{heatgg}{Dataset}(Dat, order.samples.by = NULL, facet = NULL,
      sample.id.var = "SAMPLEID", col.name = "Taxon", 
      value.name = "Abundance", trans = "log", guide = "colourbar",
      gradientn.colours = c("white","#67001F"), facet.scales = "free",
      cluster = FALSE, distfun = dist )
}

\arguments{
  \item{Tab}{A numeric matrix with samples as columns and taxons as rows}
  \item{Map}{A data.frame where the row names correspond to the column names in Tab, matching metadata to eacj sample}
  \item{Dat}{A Dataset object.}
  \item{order.samples.by}{If passes, it must be the name of a variable in Map (or in the Map attribute of the Dataset object). Samples will be ordered according to this variable.}
  \item{facet}{If passed, samples wll be arrranged in panes accroding to \code{\link{facet_grid}}. The variable names should be in the Map object.}
  \item{sample.id.var}{This is the name the sample axis will have in the heatmap.}
  \item{col.name}{The name the x-axis will have in the heatmap}
  \item{value.name}{The name that will eb on the legend indicating the meaning of the color scale}
  \item{trans}{Transformation to apply to the data during plotting. Ay continuous tranformation supported by ggplot2 can be used. Normally "trans" or "sqrt"}
  \item{guide}{Type of legend to use. This will be passed to the guide option of \code{\link{scale_fill_gradientn}}. Can be a string like "colourbar" or a function like \code{\link{guide_colourbar}}}
  \item{gradientn.colours}{Colour range to use in the heatmap.}
  \item{facet.scales}{Argument to pass to the scales option of \code{\link{facet_grid}}. Only used if option facet is not NULL.}
  \item{cluster}{Logical. Indicates wheather to perform 2 dimensional hierarchichal clustering or not.}
  \item{distfun}{Distance function to use in the hierarchichal clustering. Only relevant if cluster equeals TRUE.}
}

\details{Currently, it supports only 2d clustering. While the faceting and order options are only available for samples and when cluster is FALSE.}

\value{If cluster is FALSE it will return a ggplot2 object that can be furhter manipulated and replotted.

If cluster is TRUE it will return a list containing the heatmap and each of the dendogram plots separately
}

\references{
https://github.com/chr1swallace/random-functions/blob/master/R/ggplot-heatmap.R
http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
http://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
}

\author{Sur from Dangl Lab}

\examples{
data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)
Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)
heatgg(Tab = Dat$Tab, Map = Dat$Map, order.samples.by = "fraction")
heatgg(Tab = Dat$Tab, Map = Dat$Map, facet = fraction ~ .,
        order.samples.by = "accession")
heatgg(Tab = Dat$Tab, Map = Dat$Map, facet = fraction + accession ~ .)
heatgg(Tab = Dat$Tab, Map = Dat$Map, facet = fraction + accession ~ .,
        gradientn.colours = c("white", "black"))
p1 <- heatgg(Tab = Dat$Tab, Map = Dat$Map,
    gradientn.colours = c("white","black"),cluster = TRUE)

heatgg(Dat = Dat, order.samples.by = "fraction")
heatgg(Dat = Dat, facet = fraction ~ ., order.samples.by = "accession")
heatgg(Dat = Dat, facet = fraction + accession ~ .,)
heatgg(Dat = Dat, facet = fraction + accession ~ .,
        gradientn.colours = c("white", "black"))
p1 <- heatgg(Dat = Dat, gradientn.colours = c("white","black"),
      cluster = TRUE)
}
