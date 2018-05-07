#' Plot single variable.
#' 
#' Make a box and whiskers plot of a single variable,
#' while grouping by a variable.
#'
#' @param Map A data.frame where each row must correspond
#' to a sample, and the row names represent sample IDs.
#' @param Dat A Dataset object
#' @param var.name String indicating which variable to plot.
#' Should correspond to a header name in Map
#' @param x String indicating which variable to use as the
#' x-axis in the plot. Should correspond to a header name in Map
#' @param col String indicating which variable to use to
#' color the plot. Can be the same or different as x. It should
#' correspond to a header name in Map.
#' @param theme A ggplot2 theme to be used woth the plot.
#'
#' @return A ggplot2 plot.
#' @export
#' @author Sur Herrera Paredes
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' data(Rhizo.tax)
#' Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)
#' Dat$Map$Richness <- colSums(Dat$Tab > 0)
#' 
#' # Passing directly a data_frame
#' plotgg_var(Dat$Map, var.name = "Richness",
#'            x = "fraction")
#' plotgg_var(Dat$Map, var.name = "Richness",
#'            x = "fraction",col = "accession")
#' 
#' # Passing directly a Dataset object
#' plotgg_var(Dat, var.name = "Richness",
#'            x = "fraction")
#' plotgg_var(Dat, var.name = "Richness",
#'            x = "fraction",
#'            col = "accession")
plotgg_var <- function(...) UseMethod("plotgg_var")

#' @rdname plotgg_var
#' @method plotgg_var default
#' @export
plotgg_var.default <- function(Map, var.name ,x,
                               col = NULL,
                               theme = theme_blackbox()){
  
  p1 <- ggplot(Map, aes_string(x = x, y = var.name,
                               col = col, fill = col)) +
    geom_boxplot(fill=NA, outlier.colour = NA,
                 position = position_dodge(width = 0.9),
                 size = 2)
  
  if(is.null(col)){
    p1 <- p1 +
      geom_point(aes(fill = "black"),
                 position = position_jitterdodge(dodge.width = 0.9,
                                                 jitter.width = 0.1),
                 size = 4, shape = 21, col = "black")
  }else{
    p1 <- p1 +
      geom_point(position = position_jitterdodge(dodge.width = 0.9,
                                                 jitter.width = 0.1),
                 size = 4, shape = 21, col = "black")
  }
  
  p1 <- p1 + theme  
  return(p1)
}

#' @rdname plotgg_var
#' @method plotgg_var Dataset
#' @export
plotgg_var.Dataset <- function(Dat, var.name ,
                               x, col = NULL,
                               theme = theme_blackbox()){
  
  p1 <- plotgg_var.default(Map = Dat$Map,
                           var.name =  var.name,
                           x = x,
                           col = col,
                           theme = theme)
  
  return(p1)
}