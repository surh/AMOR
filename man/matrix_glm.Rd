\name{matrix_glm}
\alias{matrix_glm}
\alias{matrix_glm.default}
\alias{matrix_glm.Dataset}

\title{Fit Generalized Linear Models on a matrix of observations}

\description{Takes a matrix of observations, and fits the specified GLM on each species, returning the matrix of coefficients.}

\usage{
  \method{matrix_glm}{default}(x,Map,formula,family=poisson(link="log"),
                                response.name="Count",verbose=FALSE,...)
  \method{matrix_glm}{Dataset}(x,formula,family=poisson(link="log"),
                                response.name="Count",verbose=FALSE,...)
 }

\arguments{
  \item{x}{ Either a matrix of abundances with samples as columns and species as rows or a \code{Dataset} object}
  \item{Map}{A data.frame containing the variables to be modelled as columns and samples as rows. The rows should be named with sample IDs and must correspond to the column names from x if an abundance matrix was passed}
  \item{formula}{A formula specifying the model to be fit. Only the right hand side of the equation must be passed, and anything on the left side wild be silently ignored}
  \item{family}{The GLM fmaily to be used. It must be an object of class \code{family}. use \code{?family} for more help}
  \item{response.name}{String indicating the name to be used for the response (dependent) variable in the GLM}
  \item{verbose}{Logical value. If true the taza name is printed while the fit is in progress.}
  \item{...}{Other parameters for the glm() function like control.}
}

\details{This function will take each species (row) in the abundance matrix and fit the specified (cia the formula option) GLM independently on each one.}

\value{
Returns a \code{matrix.glm} object which is a list containing the following elements:
  \item{coefficients}{A matrix of coefficients for the glm fit. The matix has dimensions S x p, where S is the number of species in the abundance matrix, and p the number of parameters in the specified GLM. Each row corresponds to an independent fit of the model with glm() function and the specified formula}
  \item{SE}{Similar to Coef, but this matrix contains the estimated Standard Errors (ie. the Standard deviation of the coefficients).}
  \item{AIC}{A named vector containing the Akaike Information Criteria (AIC) of each independent fit. The names of the vector correspond to the species IDs in the input abundance matrix}
  \item{call}{Function call via match.call}
  \item{family}{GLM model family. An object of class \code{family}}
  \item{X}{Design matrix of the model fit.}
}

\author{
Sur from Dangl Lab
}

