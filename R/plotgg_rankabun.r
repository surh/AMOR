plotgg_rankabun <- function(...) UseMethod("plotgg_rankabun")

plotgg_rankabun.default <- function(Tab, Map = NULL, groupby = NULL,
                                    sortby = NULL,
                                    theme = theme_blackbox,
                                    variable.name = "Taxon",
                                    value.name = "Abundance",
                                    FUN = mean){
#   Tab <- Dat$Tab
#   Map <- Dat$Map
#   groupby <- "fraction"
#   sortby <- "E"
#   FUN <- mean
#   variable.name <- "Taxon"
#   value.name <- "Abundance"
#   theme <- theme_blackbox
  
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

plotgg_rankabun.Dataset <- function(Dat, groupby = NULL,
                                    sortby = NULL,
                                    theme = theme_blackbox,
                                    variable.name = "Taxon",
                                    value.name = "Abundance",
                                    FUN = mean){
  res <- plotgg_rankabun.default(Tab = Dat$Tab, Map = Dat$Map, groupby = groupby, sortby = sortby, theme = theme,
                          variable.name = variable.name, value.name = value.name, FUN = FUN)
  return(res)
}