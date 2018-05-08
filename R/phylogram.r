#' Phylogram
#' 
#' Make a relative abundance barplot of taxon abundances.
#'
#' @param Tab A matrix object representing samples as columns and taxa as rows.
#' @param Dat A Dataset object.
#' @param Map A data frame with metadata for Tab. One row per sample, and the
#' rows must be named to match the column names in Tab, and they must be in the same order.
#' @param facet facet formula for facet_grid.
#' @param colname Label to be used in the x-axis of the phylogram.
#' @param variable.name Label to be used in for the colors in the phylogram
#' @param value.name Label for the y-axis of the phylogram.
#' @param scales scales option for \link{facet_grid}.
#' @param space space option for \link{facet_grid}
#' @param nrow.legend number of rows to use in the legend.
#' @param ntaxa Number of taxa to plot, the rest will be collapsed into
#' one category. Samples will be sorted by abundance and the top ntaxa will be plotted.
#' @param other_name Name to give to the collapsed taxa when ntaxa is used.
#'
#' @return A ggplot2 object of the plot.
#' @export
#' @author Sur Herrera Paredes
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' data(Rhizo.tax)
#' Dat <- create_dataset(Rhizo,Rhizo.map,Rhizo.tax)
#' 
#' Dat.phyl <- collapse_by_taxonomy(Dat = Dat, level = 4)
#' phylogram(Tab = Dat.phyl$Tab,
#'           Map = Dat.phyl$Map,
#'           facet = ~ fraction)
#' phylogram(Dat.phyl, facet = ~ fraction)
#' phylogram(Dat.phyl, facet = ~ fraction,
#'           ntaxa = 5)
phylogram <- function(...) UseMethod("phylogram")

#' @rdname phylogram
#' @method phylogram default
#' @export
phylogram.default <- function(Tab,Map=NULL,facet = NULL,colname="Sample",
                              variable.name="Taxon",value.name="Abundance",
                              scales="free_x",space = "free_x", nrow.legend=20,
                              ntaxa = NULL, other_name = "other"){
  
  # Collapse extra taxa if needed
  if(is.numeric(ntaxa)){
    if(nrow(Tab) > ntaxa){
      # Sort by mean abundance
      select <- rowSums(Tab)
      select <- sort(select,decreasing = TRUE)
      select <- names(select)[1:ntaxa]
      groups <- row.names(Tab)
      groups[ !(groups %in% select) ] <- other_name
      Tab <- collapse_matrix(x = Tab,groups = groups,dim = 1,FUN = sum)
    }
  }
  
  # Add Sample id column
  Tab <- as.data.frame(t(Tab))
  measure.vars <- names(Tab)
  Tab[,colname] <- row.names(Tab)
  
  # If needed and provided add mapping file for faceting
  if(!is.null(Map) & !is.null(facet)){
    Map <- Map[row.names(Tab),]
    colname <- c(colname,names(Map))
    Tab <- cbind(Tab,Map)
  }
  
  # Prepare data for plotting
  Dat <- melt(Tab,id.vars=colname,
              variable.name=variable.name,
              value.name=value.name)
  
  # Plot
  p1 <- ggplot(Dat,aes_string(x=colname,y=value.name,fill=variable.name)) +
    geom_bar(stat = "identity",position = "fill", width = 1) +
    theme(axis.text = element_text(angle=90))
  
  if(!is.null(Map) & !is.null(facet)){
    p1 <- p1 + facet_grid(facet,scales = scales,space = space)
  }
  
  # Modify legend
  p1 <- p1 + guides(fill=guide_legend(nrow=nrow.legend))
  
  return(p1)
}

#' @rdname phylogram
#' @method phylogram Dataset
#' @export
phylogram.Dataset <- function(Dat, facet = NULL,colname="Sample",
                              variable.name="Taxon",value.name="Abundance",
                              scales="free_x",space = "free_x", nrow.legend=20,
                              ntaxa = NULL,other_name = "other"){
  
  res <- phylogram.default(Tab = Dat$Tab, Map = Dat$Map,
                           facet = facet, colname = colname,
                           variable.name = variable.name,
                           value.name = value.name, scales = scales,
                           space = space, nrow.legend = nrow.legend,
                           ntaxa = ntaxa,
                           other_name = other_name)
  
  return(res)
}