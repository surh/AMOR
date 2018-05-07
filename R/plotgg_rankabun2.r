#' Make rank abundance plot with Poisson confidence intervals.
#' 
#' Make a rank abundance plot with Poisson confidence intervals
#' from a Dataset object.
#' IMPORTANT: There is no guarantee that Poisson confidence
#' intervals are representative of the data. This should be viewed
#' only as descriptive statistic and not for inference.
#' 
#' 
#' 
#' 
#' @param Tab A matrix object with samples as columns and taxa as rows.
#' @param Map A data frame with metadata forTab. Each row must be a
#' sample with row names matching column names in Tab and in the
#' same order as in Tab.
#' @param Dat A Dataset object.
#' @param groupby Variable name to be used for grouping samples before plotting th rank abundance.
#' @param sortby Variable value to be used as reference for sorting the taxa in the rank abundance plot.
#' @param alpha Transparency parameter for ggplot2
#' @param theme ggplot2 theme to be used for plotting.
#' @param variable.name x-axis label in the plot.
#' @param value.name y-axis label in the plot.
#' @param sample.id.name name to store sample IDs in Map. Used for internal handling only.
#'
#' @return A ggplot2 plot.
#' @export
#' @author Sur Herrera Paredes
#' @seealso \link{plotgg_rankabun}
#'
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' data(Rhizo.tax)
#' Dat <- create_dataset(Rhizo, Rhizo.map, Rhizo.tax)
#' 
#' plotgg_rankabun2(Tab = Dat$Tab,Map = Dat$Map,
#'                  groupby = "fraction",sortby = "E")
#' plotgg_rankabun2(Dat, groupby = "fraction", sortby = "E")
plotgg_rankabun2 <- function(...) UseMethod("plotgg_rankabun2")

#' @rdname plotgg_rankabun2
#' @method plotgg_rankabun2 default
plotgg_rankabun2.default <- function(Tab, Map, groupby, sortby,
                                     alpha = 0.2, theme = theme_blackbox,
                                     variable.name = "Taxon",
                                     value.name = "Abundance",
                                     sample.id.name = "SAMPLEID"){
  # This is an experimental function. SHOULD NOT BE USED AS A FORMAL TEST OF SIGNIFICANCE.
  
  if(any(colnames(Tab) != row.names(Map))){
    stop("ERROR: names in Tab and Map do not match",call. = TRUE)
  }
  
  original_variables <- colnames(Map)
  if(sample.id.name %in% original_variables){
    stop("ERROR: The sample.id.name must not be already present as a variable in Map",call. = TRUE)
  }
  Dat <- cbind(Map,t(Tab))
  Dat[,sample.id.name] <- row.names(Dat)
  
  # Convert data to ggplot format
  Dat <- melt(Dat,id.vars = c(original_variables,sample.id.name),
              variable.name = variable.name,value.name = value.name)
  
  # Group and calculate poisson CIs
  Dat2 <- aggregate(Dat[,value.name],by = list(Dat[,variable.name],Dat[,groupby]),sum)
  colnames(Dat2) <- c(variable.name, groupby, value.name)
  Nsamples <- table(Map[,groupby])
  Dat2$N <- Nsamples[ as.character(Dat2[,groupby]) ]
  Dat2$upper <- Dat2[, value.name] + 1.96*sqrt(Dat2[ , value.name ])
  Dat2$lower <- Dat2[ , value.name ] - 1.96*sqrt(Dat2[ , value.name ])
  Dat2[,value.name] <- Dat2[, value.name] / Dat2$N
  Dat2$lower <- Dat2$lower / Dat2$N
  Dat2$upper <- Dat2$upper / Dat2$N
  
  # Order
  order_vector <- subset(Dat2, Dat2[ , groupby ] == sortby)
  order_vector <- order_vector[ order(order_vector[, value.name ],decreasing = TRUE), ]
  order_vector <- unique(as.character(order_vector[, variable.name]))
  Dat2[, variable.name] <- factor(Dat2[, variable.name], levels = order_vector)
 
  # Format
  Dat2[,value.name] <- as.numeric(Dat2[,value.name])
  Dat2[,"lower"] <- as.numeric(Dat2[,"lower"])
  Dat2[,"upper"] <- as.numeric(Dat2[,"upper"])
  
  p1 <- ggplot(Dat2,aes(x = Taxon,
                        y = Abundance,
                        group=fraction,
                        col = fraction,
                        fill=fraction)) +
    geom_line() +
    geom_ribbon(aes(x = as.numeric(Dat2$Taxon),
                    ymin=lower,ymax=upper), alpha = 0.2)
  p1
  
  # Plot
  p1 <- ggplot(Dat2, aes_string(x = variable.name,
                                y = value.name,
                                col = groupby,
                                group = groupby,
                                fill = groupby))+
    geom_line() +
    geom_ribbon(aes(x=as.numeric(Dat2$Taxon),ymin=lower,ymax=upper), alpha=0.2) +
    theme
  p1
  

  
  return(p1)
}

#' @rdname plotgg_rankabun2
#' @method plotgg_rankabun2 Dataset
plotgg_rankabun2.Dataset <- function(Dat, groupby, sortby,
                                     alpha = 0.2, theme = theme_blackbox,
                                     variable.name = "Taxon",
                                     value.name = "Abundance",
                                     sample.id.name = "SAMPLEID"){
  
  p1 <- plotgg_rankabun2.default(Tab = Dat$Tab,
                                 Map = Dat$Map,
                                 groupby = groupby,
                                 sortby = sortby,
                                 alpha = alpha ,
                                 theme = theme,
                                 variable.name = variable.name,
                                 value.name = value.name,
                                 sample.id.name = sample.id.name)
  return(p1)
}
