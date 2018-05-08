#' Plot single taxon.
#' 
#' Make a box and whiskers plot of a single taxon, while grouping by a variable.
#'
#' @param Tab Matrix with samples as columns and taxa as rows.
#' @param Dat A Dataset object.
#' @param Map Data frame with metadata for Tab. Each row must
#' correspond to a sample in Tab, they must be in the same order,
#' and be named with the same names as the columns in Tab
#' @param taxon String indicating which taxon to plot.
#' @param x String indicating which variable to use as the x-axis
#' in the plot. Should correspond to a header name in Map
#' @param col String indicating which variable to use to color the
#' plot. Can be the same or different as x. It should correspond to
#' a header name in Map.
#' @param var.name y-axis label on the plot.
#' @param theme A ggplot2 theme to be used woth the plot.
#'
#' @return A ggplot2 plot
#' @export
#' @author Sur Herrera Paredes
#' @seealso \link{plotgg_var}
#' 
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' data(Rhizo.tax)
#' Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)
#' 
#' # Passing independently an abundance table and mapping file
#' plotgg_taxon(Tab = Dat$Tab, Map = Dat$Map,
#'              taxon = "OTU_14834", x = "fraction")
#' plotgg_taxon(Tab = Dat$Tab, Map = Dat$Map,
#'              taxon = "OTU_14834", x = "fraction",
#'              col = "accession")
#' 
#' # Passing a dataset object
#' plotgg_taxon(Dat, taxon = "OTU_14834",
#'              x = "fraction")
#' plotgg_taxon(Dat, taxon = "OTU_14834",
#'              x = "fraction",
#'              col = "accession")
plotgg_taxon <- function(...) UseMethod("plotgg_taxon")

#' @rdname plotgg_taxon
#' @method plotgg_taxon default
#' @export
plotgg_taxon.default <- function(Tab, Map,
                                 taxon ,x,
                                 col = NULL,
                                 var.name = "Abundance",
                                 theme = theme_blackbox()){
  
  if(var.name %in% names(Map)){
    stop("ERROR: var.name exists already in Map",call. = TRUE)
  }
  
  if(all(colnames(Tab) == row.names(Map))){
    Dat <- Map
    Dat[,var.name] <- as.vector(Tab[ taxon, ])
    
    p1 <- ggplot(Dat,aes_string(x = x, y = var.name,
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
  }else{
    stop ("ERROR: Sample names in Tab do not match sample names in Map\n",
          call. = TRUE)
  }
  
  p1 <- p1 + theme
  
  return(p1)
}

#' @rdname plotgg_taxon
#' @method plotgg_taxon Dataset
#' @export
plotgg_taxon.Dataset <- function(Dat, taxon , x,
                                 col = NULL,
                                 var.name = "Abundance",
                                 theme = theme_blackbox()){
  
  p1 <- plotgg_taxon.default(Tab = Dat$Tab, Map = Dat$Map,
                             taxon = taxon , x = x, col = col,
                             var.name = var.name, theme = theme)
  
  return(p1)
}
