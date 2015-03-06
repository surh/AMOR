\name{create_dataset}
\alias{create_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Dataset object
}
\description{
Takes an abundance matrix, a metadata data.frame and (optionally) a taxonomy definition, and creates a Dataset object.
}
\usage{
create_dataset(Tab=NULL,Map=NULL,Tax=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tab}{A numeric matrix (or a numerci object coercible into a matrix) of S x N dimensions. Whenre S is the number of different taxa (OTUs, taxonomic levels etc), and N the number of samples. The rows and columns must be named with the ID of taxons and samples respectively}
  \item{Map}{Optional data.frame of dimensions N x p, where N is the number of samples and p the number of variables for which there is information. The rows must be named with the sample IDs, and these IDs must correspond to the same order as in the Tab argument}
  \item{Tax}{Optional data.frame of dimensions S x 2, where S is the number of taxons. The rows. must be named according to the taxon IDs, which should correspond to the txons in the Tab argument. The dataframe must contain columns ID and Taxonomy, which contain the Taxa ID and the taxonomy string, in the RDP classifier format}
}
\details{See description}
\value{
Returns a list of class \code{Dataset} with the following elements:
  \item{Tab}{ A numeric matrix containing counts of species (rows) per sample (columns) }
  \item{Map}{ A data.frame containing the Metadata for samples in Tab }
  \item{Tax}{ Optional data.frame containing taxonomic annotation for species in Tab }
}
\author{Sur from Dangl Lab}
