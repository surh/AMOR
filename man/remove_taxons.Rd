\name{remove_taxons}
\alias{remove_taxons}

\title{Remove taxa from a Dataset object.}

\description{Remove taxa from a Dataset object.}

\usage{
remove_samples(Dat,taxons)
}

\arguments{
  \item{Dat}{A Dataset object}
  \item{taxons}{String vector indicating which taxa to remove. Taxon IDs must be equal to the row names of the abundance table on \code{Dat}.}
}

\details{}

\value{A Dataset object.}

\author{Sur from Dangl Lab}

\seealso{\link{remove_samples}}

\examples{
data(Rhizo)
data(Rhizo.map)
Dat <- create_dataset(Rhizo,Rhizo.map)
Dat <- remove_samples(Dat,colnames(Dat$Tab)[1:50])
Dat <- remove_taxons(Dat,row.names(Dat$Tab)[1:115])
sum(colSums(Dat$Tab) == 0)
sum(rowSums(Dat$Tab) == 0)
Dat$Tab
Dat.clean <- clean(Dat)
Dat.clean$Tab
}