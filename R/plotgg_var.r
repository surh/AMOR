plotgg_var <- function(...) UseMethod("plotgg_var")

plotgg_var.default <- function(Map, var.name ,x, col = NULL, theme = theme_blackbox){
  # Map <- Map
  # x <- "Fraction"
  # col <- "Genotype"
  # var.name <- "Depth"
  
  p1 <- ggplot(Map,aes_string(x = x, y = var.name, col = col, fill = col)) +
    geom_boxplot(fill=NA,outlier.colour = NA, position = position_dodge(width = 0.9), size = 2)
  
  if(is.null(col)){
    p1 <- p1 +
      geom_point(aes(fill = "black"),position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), size = 4, shape = 21, col = "black")
  }else{
    p1 <- p1 +
      geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), size = 4, shape = 21, col = "black")
  }
  
  p1 <- p1 + theme  
  return(p1)
}


plotgg_var.Dataset <- function(Dat, var.name ,x, col = NULL, theme = theme_blackbox){
  p1 <- plotgg_var.default(Map = Dat$Map, var.name =  var.name, x = x, col = col, theme = theme)
  
  return(p1)
}
