\name{plotgg.PCO}
\alias{plotgg.PCO}

\title{Plotting a PCoA}

\description{Function for plotting results of \code{\link{PCO}}.}

\usage{
  \method{plotgg}{PCO}(x,components=c("V1","V2"),shape=NULL,col=NULL,
                       point_size=2)
}

\arguments{
  \item{x}{A \code{PCA} object.}
  \item{components}{Vector of length 2 indicating which components to plot.}
  \item{shape}{String indicating which variable to use as aestetics mapping for shape. Must correspond to a column header in the Map attribute of the PCO object.}
  \item{col}{String indicating which variable to use as aestetics mapping for color. Must correspond to a column header in the Map attribute of the PCO object.}
  \item{point_size}{size for the points in the plot}
}

\details{}

\value{A ggplot2 object of the PCoA plot}

\author{Sur from Dangl Lab}

\seealso{\link{PCO}}

\examples{
  data(Rhizo)
  data(Rhizo.map)
  Dat <- create_dataset(Rhizo,Rhizo.map)
  
  # distfun <- function(x) vegan::vegdist(x,method="bray") #requires vegan package
  distfun <- dist
  
  Dat.pco <- PCO(Dat,dim=2,distfun=distfun)
  plotgg(Dat.pco)
  plotgg(Dat.pco,shape="fraction",point_size=3)
  plotgg(Dat.pco,shape="fraction",col="accession",point_size=4)
}