#' Collapse matrix
#' 
#' Takes a two-dimensional matrix and collapses row or columns
#' according to a specified function.
#' 
#' This function uses plyr functions to aggregate columnwise or
#' rowwise according to a grouping factor.
#' 
#' @param x A \code{matrix} object.
#' @param groups A vector or factor of the same length as the rows or
#' columns of x, which indicates which elements to collapse.
#' @param dim Which dimension to collapse. 1 for rows and 2 for columns
#' @param FUN Function (or function name) to apply to elements per groups
#' 
#' @return Returns a \code{matrix} object.
#' 
#' @author Sur Herrera Paredes
#' 
#' @export
#' 
#' @importFrom plyr aaply
#' @importFrom plyr daply
#' 
#' @examples 
#' x <- cbind(rbind(matrix(0,2,2),matrix(1,2,2)),
#'            rbind(matrix(1,2,2),matrix(0,2,2)))
#' colnames(x) <- c("C1","C2","C3","C4")
#' row.names(x) <- c("R1","R2","R3","R4")
#' 
#' x
#' collapse_matrix(x=x,groups = rep(c(1,2),each=2),dim=1,FUN=sum)
#' collapse_matrix(x=x,groups = rep(c(1,2),each=2),dim=2,FUN=mean)
collapse_matrix <- function(x,groups,dim=1,FUN=sum){

  if(!is.matrix(x)){
    # From documentation, this checks that object is 2d matrix.
    # If ever there is extension to more dims, need to change, maybe with inherits or isa.
    stop("ERROR: x must be a matrix ",call.=TRUE)
  }
  if(length(groups) != length(dimnames(x)[[dim]])){
    stop("ERROR: length of groups does not correspond to number of vectors",call.=TRUE)
  }
  
  FUN <- match.fun(FUN)
  groups <- factor(groups)  
  
  # Tab.collapsed <- apply(x,((dim %% 2) + 1),tapply,
  #                        INDEX = groups, FUN)
  Tab.collapsed <- plyr::aaply(.data = x, .margins = ((dim %% 2) + 1),
                               .fun = function(x, g){
                                 d <- data.frame(x = x, g = g)
                                 plyr::daply(.data = d, .variables = ~ g,
                                             .fun = function(x){FUN(x$x)})
                               }, g = groups, .drop = FALSE)
  
  if(dim == 1){
    Tab.collapsed <- t(Tab.collapsed)
  }
  
  return(Tab.collapsed)
}
