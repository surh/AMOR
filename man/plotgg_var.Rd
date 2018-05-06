\name{plotgg_var}
\alias{plotgg_var}
\alias{plotgg_var.default}
\alias{plotgg_var.Dataset}

\title{Plot single variable.}

\description{Make a box and whiskers plot of a single variable, while grouping by a variable.}

\usage{
  \method{plotgg_var}{default}(Map, var.name, x, col = NULL,
                        theme = theme_blackbox)
  \method{plotgg_var}{Dataset}(Dat, var.name, x, col = NULL,
                        theme = theme_blackbox)
}

\arguments{
  \item{Map}{A data.frame where each row must correspond to a sample, and the row names represent sample IDs.}
  \item{Dat}{A Dataset object.}
  \item{var.name}{String indicating which variable to plot. Should correspond to a header name in \code{Map}}
  \item{x}{String indicating which variable to use as the x-axis in the plot. Should correspond to a header name in \code{Map}}
  \item{col}{String indicating which variable to use to color the plot. Can be the same or different as \code{x}. It should correspond to a header name in \code{Map}.}
  \item{theme}{A ggplot2 theme to be used woth the plot.}
}

\value{A ggplot2 plot.}

\author{Sur from Dangl Lab}

\examples{
data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)
Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)
Dat$Map$Richness <- colSums(Dat$Tab > 0)

# Passing directly a data_frame
plotgg_var(Dat$Map,var.name = "Richness",x = "fraction")
plotgg_var(Dat$Map,var.name = "Richness",x = "fraction",col = "accession")

# Passing directly a Dataset object
plotgg_var(Dat,var.name = "Richness",x = "fraction")
plotgg_var(Dat,var.name = "Richness",x = "fraction",col = "accession")
}
