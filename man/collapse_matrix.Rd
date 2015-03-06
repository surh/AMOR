\name{collapse_matrix}
\alias{collapse_matrix}

\title{Collapse matrix}

\description{Takes a two-dimensional matrix and collapses row or columns according to a specified function}

\usage{collapse_matrix(x,groups,dim=1,FUN=sum)}

\arguments{
  \item{x}{A \code{matrix} object.}
  \item{groups}{A vector or factor of the same length as the rows or columns of x, which indicates which elements to collapse}
  \item{dim}{Which dimension to collapse. 1 for rows and 2 for columns}
  \item{FUN}{Function (or function name) to apply to elements per groups.}
}

\details{}

\value{Returns a \code{matrix} object.}

\author{Sur from the Dangl Lab.}

\examples{
x <- cbind(rbind(matrix(0,2,2),matrix(1,2,2)),
            rbind(matrix(1,2,2),matrix(0,2,2)))
colnames(x) <- c("C1","C2","C3","C4")
row.names(x) <- c("R1","R2","R3","R4")

x
collapse_matrix(x=x,groups = rep(c(1,2),each=2),dim=1,FUN=sum)
collapse_matrix(x=x,groups = rep(c(1,2),each=2),dim=2,FUN=mean)
}
