plotgg_taxon <- function(...) UseMethod("plotgg_taxon")

plotgg_taxon.default <- function(Tab, Map,
                                 taxon ,x,
                                 col = NULL,
                                 var.name = "Abundance",
                                 theme = theme_blackbox){
  
  if(var.name %in% names(Map)){
    stop("ERROR: var.name exists already in Map",call. = TRUE)
  }
  
  if(all(colnames(Tab) == row.names(Map))){
    Dat <- Map
    Dat[,var.name] <- as.vector(Tab[ taxon, ])
    
    p1 <- ggplot(Dat,aes_string(x = x, y = var.name, col = col, fill = col)) +
      geom_boxplot(fill=NA,outlier.colour = NA, position = position_dodge(width = 0.9), size = 2)
    
    if(is.null(col)){
      p1 <- p1 +
        geom_point(aes(fill = "black"),position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), size = 4, shape = 21, col = "black")
    }else{
      p1 <- p1 +
        geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), size = 4, shape = 21, col = "black")
    }
  }else{
    stop ("ERROR: Sample names in Tab do not match sample names in Map\n",call. = TRUE)
  }
  
  p1 <- p1 + theme
  
  return(p1)
}

plotgg_taxon.Dataset <- function(Dat, taxon , x,
                                 col = NULL,
                                 var.name = "Abundance",
                                 theme = theme_blackbox){
  p1 <- plotgg_taxon.default(Tab = Dat$Tab, Map = Dat$Map, taxon = taxon , x = x, col = col, var.name = var.name,theme = theme)
  return(p1)
}

