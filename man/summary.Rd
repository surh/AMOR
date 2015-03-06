\name{summary}
\alias{summary.matrix.glm}
\alias{summary.bootglm}
\alias{print.summary.matrix.glm}
\alias{print.summary.bootglm}

\title{Summarize matrix GLM results}

\description{Takes a \code{matrix.glm} or \code{bootglm} object and reformats the results}

\usage{
  \method{summary}{matrix.glm}(object,sortby="Variable",...)
  \method{summary}{bootglm}(object,sortby="Variable",...)
  
  \method{print}{summary.matrix.glm}(x,...)
  \method{print}{summary.bootglm}(x,...)
}

\arguments{
  \item{object}{A \code{matrix.glm} or \code{bootglm} object.}
  \item{sortby}{A string indicating how to sort the coefficients. Either by "Vatiable" or "Taxon"}
  \item{x}{A \code{summary.matrix.glm} or \code{summary.bootglm} object to be printed}
  \item{...}{Further arguments to summary and print methods. NOTE: ignored at this point}
}

\value{
Returns a \code{summary.matrix.glm} or \code{summary.bootglm} object which is a list containing the following elements:
  \item{coefficients}{Coefficient table. For matrix.glm objects it includes standard errors, z-values and p-values. NOTE: Since only p-values from z-tests are provided at this point, they are only appropriate for binomial, poisson and negativee.binomial models. For boot.glm objects in includes the 2.5\% and 97.5\% percentiles and the bootstrap confidence for each parameter (Pr(param > 0)).}
  \item{call}{Original matrix_glm() or bootstrap_glm() call}
}

\author{
Sur from Dangl Lab
}

