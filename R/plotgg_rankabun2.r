plotgg_rankabun2 <- function(...) UseMethod("plotgg_rankabun2")

plotgg_rankabun2.default <- function(Tab, Map, groupby, sortby, alpha = 0.2, theme = theme_blackbox,
                                    variable.name = "Taxon", value.name = "Abundance", sample.id.name = "SAMPLEID"){
  # This is an experimental function. SHOULD NOT BE USED AS A FORMAL TEST OF SIGNIFICANCE.
#   Tab <- Dat$Tab
#   Map <- Dat$Map
#   variable.name <- "Taxon"
#   value.name <- "Abundance"
#   sample.id.name <- "SAMPLEID"
#   alpha <- 0.2
#   groupby <- "fraction"
#   sortby <- "Soil"
#   theme <- theme_blackbox
  
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
  Dat <- melt(Dat,id.vars = c(original_variables,sample.id.name),variable.name = variable.name,value.name = value.name)
  
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
  
  # Plot
  p1 <- ggplot(Dat2, aes_string(x = variable.name, y = value.name, col = groupby, group = groupby, fill = groupby))+
    geom_line() +
    geom_ribbon(aes(ymin=lower,ymax=upper),alpha=alpha) +
    theme
  
  return(p1)
}

plotgg_rankabun2.Dataset <- function(Dat, groupby, sortby, alpha = 0.2, theme = theme_blackbox,
                                     variable.name = "Taxon", value.name = "Abundance", sample.id.name = "SAMPLEID"){
  p1 <- plotgg_rankabun2.default(Tab = Dat$Tab,Map = Dat$Map, groupby = groupby, sortby = sortby, alpha = alpha , theme = theme,
                                 variable.name = variable.name, value.name = value.name, sample.id.name = sample.id.name)
  return(p1)
}
