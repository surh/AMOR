#' Make rank abundance plot
#' 
#' Make a rank abundance plot from a Dataset object.
#'
#' @param Tab A matrix object with samples as columns and taxa as rows.
#' @param Map A data frame with metadata for Tab. Each row must be a
#' sample with row names matching column names in Tab and in the same order as in Tab.
#' @param Dat A Dataset object.
#' @param groupby Variable name to be used for grouping samples before plotting th rank abundance.
#' @param sortby Variable value to be used as reference for sorting the taxa in the rank abundance plot.
#' @param theme ggplot2 theme to be used for plotting.
#' @param variable.name x-axis label in the plot.
#' @param value.name y-axis label in the plot.
#' @param FUN Function to use to aggregate samples according to groupby
#'
#' @return A ggplot2 plot.
#' @export
#' @author Sur Herrera Paredes
#' @seealso \link{plotgg_rankabun2}
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' data(Rhizo.tax)
#' Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)
#' 
#' plotgg_rankabun(Tab = Dat$Tab)
#' plotgg_rankabun(Tab = Dat$Tab, Map = Dat$Map,groupby = "fraction",sortby = "Soil")
#' 
#' plotgg_rankabun(Dat)
#' plotgg_rankabun(Dat, groupby = "fraction", sortby = "E")
plotgg_rankabun <- function(...) UseMethod("plotgg_rankabun")

#' @rdname plotgg_rankabun
#' @method plotgg_rankabun default
#' @export
plotgg_rankabun.default <- function(Tab, Map = NULL, groupby = NULL,
                                    sortby = NULL,
                                    theme = theme_blackbox(),
                                    variable.name = "Taxon",
                                    value.name = "Abundance",
                                    FUN = mean){
  
  Dat <- create_dataset(Tab,Map)
  
  if(!is.null(groupby) && !is.null(sortby)){
    Dat <- pool_samples.Dataset(x = Dat,groups = groupby,FUN = FUN, return.dataset = TRUE)
    order_vector <- row.names(Dat$Tab)[ order(Dat$Tab[,sortby], decreasing = TRUE) ]
  }else{
    order_vector <- apply(Dat$Tab,1,FUN)
    order_vector <- sort(order_vector,decreasing = TRUE)
    order_vector <- names(order_vector)
    
  }
  
  Dat2 <- data.frame(Dat$Tab)
  Dat2[,variable.name] <- row.names(Dat$Tab)
  Dat2 <- melt(Dat2,id.vars = variable.name,value.name = value.name,variable.name = "SAMPLEID")
  Dat2[,variable.name] <- factor(Dat2[,variable.name], levels = order_vector)
  
  p1 <- ggplot(Dat2,aes_string(x = variable.name, y = value.name))
  if(!is.null(groupby) && !is.null(sortby)){
    p1 <- p1 + geom_line(aes_string(group = "SAMPLEID",col = "SAMPLEID"))
  }else{
    p1 <- p1 + geom_line(aes_string(group = "SAMPLEID"))
  }
  p1 <- p1 + theme
  # p1
  return(p1)
}

#' @rdname plotgg_rankabun
#' @method plotgg_rankabun Dataset
#' @export
plotgg_rankabun.Dataset <- function(Dat, groupby = NULL,
                                    sortby = NULL,
                                    theme = theme_blackbox(),
                                    variable.name = "Taxon",
                                    value.name = "Abundance",
                                    FUN = mean){
  
  res <- plotgg_rankabun.default(Tab = Dat$Tab,
                                 Map = Dat$Map,
                                 groupby = groupby,
                                 sortby = sortby,
                                 theme = theme,
                                 variable.name = variable.name,
                                 value.name = value.name,
                                 FUN = FUN)
  return(res)
}