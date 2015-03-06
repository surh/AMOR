\name{bootstrap_glm}
\alias{bootstrap_glm}
\alias{bootstrap_glm.default}
\alias{bootstrap_glm.Dataset}

\title{Bootstrap estimation of Generalized Linear Models on a matrix of observations}

\description{Takes a matrix of observations, and fits the specified GLM on each species. It perfomrs N bootstrap pesudoreplicates and returns the matrix of coefficients, for the original dataset as well as bootstrap estimates. It uses \code{glm_matrix()} to fit each of the pseudoreplicates}

\usage{
  \method{bootstrap_glm}{default}(x,Map,formula,family=poisson(link="log"),
  	response.name="Count",N=100,verbose=FALSE)
  \method{bootstrap_glm}{Dataset}(x,formula,family=poisson(link="log"),
  	response.name="Count",N=100,verbose=FALSE)
 }

\arguments{
  \item{x}{ Either a matrix of abundances with samples as columns and species as rows or a \code{Dataset} object.}
  \item{Map}{A data.frame containing the variables to be modelled as columns and samples as rows. The rows should be named with sample IDs and must correspond to the column names from x if an abundance matrix was passed.}
  \item{formula}{A formula specifying the model to be fit. Only the right hand side of the equation must be passed, and anything on the left side wild be silently ignored.}
  \item{family}{The GLM fmaily to be used. It must be an object of class \code{family}. use \code{?family} for more help.}
  \item{response.name}{String indicating the name to be used for the response (dependent) variable in the GLM.}
  \item{N}{Number of bootstrap pseudoreplicates to perform.}
  \item{verbose}{Logical parameter. IF true each call to matrix_glm will print progress information.}
  \item{...}{Other parameters for the glm() function like control.}
}

\description{This works by calling \code{matrix_glm()} on the original data, and each pf the N bootstrap pseudoreplicates generated.}

\value{
Returns a \code{bootglm} object which is a list containing the following elements:
  \item{orig}{A \code{matrix.glm} object containing the output of \code{matrix.glm()} on the original data.}
  \item{boot}{A list containing elements: 1) coefficients, a three-dimensioanl array containing the matrix of coefficients for each bootstrap pseudoreplicate; and 2) SE, a three-dimensional array containing the matrix of standard errors for each bootstrap pseudoreplicate.}
  \item{coefficients}{A matrix of coefficients calculated as the mean accross each bootstrap pseudoreplicates.}
  \item{call}{Function call via match.call().}
  \item{family}{GLM model family. An object of class \code{family}.}
  \item{pseudoreplicates}{Matrix of dimensions (number iof samples) x (number of pseudoreplicates), where the i-column gives the sample ID, of the samples used for the i-th pseudoreplicate.}
  \item{N}{Number of bootstrap pseudoreplicates performed.}
}

\author{
Sur from Dangl Lab
}

