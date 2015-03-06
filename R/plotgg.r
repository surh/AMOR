plotgg <- function(...) UseMethod("plotgg")

plotgg.PCA <- function(x,components=c("PC1","PC2"),shape=NULL,col=NULL,
                       biplot=FALSE,biplot_color="grey21",point_size=2){
#   x <- Dat.pca
#   components <- c("PC1","PC2")
#   col <- NULL
#   shape <- NULL
#   col <- "accession"
#   shape <- "fraction"
#   biplot <- TRUE
#   biplot_color <- "grey21"
#   point_size <- 2
  
  Dat <- as.data.frame(x$scores)
  Dat$biplot <- "scores"
  Dat$ID <- row.names(Dat)
  if(biplot){
    Dat2 <- as.data.frame(x$loadings)
    Dat2$biplot <- "loadings"
    Dat2$ID <- row.names(Dat2)
    biplot_ratio <- max(abs(Dat[,components[1]])) / max(abs(Dat2[,components[1]]))
    Dat2[,components] <- Dat2[,components] * biplot_ratio
    Dat <- rbind(Dat,Dat2)
  }
  
  if(!is.null(x$Map)){
    loading_names <- row.names(Dat)[ !(row.names(Dat) %in% row.names(x$Map)) ]
    map <- data.frame(matrix(NA,nrow=length(loading_names),ncol=ncol(x$Map)))
    row.names(map) <- loading_names
    names(map) <- names(x$Map)
    map <- rbind(x$Map,map)

    Dat <- cbind(Dat,map)    
  }
  
  
  p1 <- ggplot(Dat,aes_string(x=components[1],y=components[2]))
  if(biplot){
    require(grid)
    p1 <- p1 + geom_segment(data = subset(Dat, biplot == "loadings"),aes_string(x=0,y=0,xend=components[1],yend=components[2]),
                            arrow=arrow(length = unit(0.5, "cm")),col=biplot_color)
  }
  p1 <- p1 + geom_point(data = subset(Dat, biplot == "scores"),aes_string(shape=shape,col=col),size=point_size)
  
  p1 <- p1 + theme(axis.text = element_text(color="black"),
                   axis.title = element_text(face="bold"),
                   panel.background = element_rect(color="black",size=3,fill=NA),
                   panel.grid = element_blank())
  #p1
    
  return(p1) 
}


plotgg.PCO <- function(x,components=c("V1","V2"),shape=NULL,col=NULL,
                       point_size=2){
#   x <- Dat.pco
#   components <- c("V1","V2")
#   col <- NULL
#   shape <- NULL
#   col <- "accession"
#   shape <- "fraction"
#   point_size <- 2
  
  Dat <- as.data.frame(x$points)
  Dat$ID <- row.names(Dat)
  
  if(!is.null(x$Map)){
    Dat <- cbind(Dat,x$Map)    
  }
  
  
  p1 <- ggplot(Dat,aes_string(x=components[1],y=components[2]))
  
  p1 <- p1 + geom_point(aes_string(shape=shape,col=col),size=point_size)
  
  p1 <- p1 + theme(axis.text = element_text(color="black"),
                   axis.title = element_text(face="bold"),
                   panel.background = element_rect(color="black",size=3,fill=NA),
                   panel.grid = element_blank())
  #p1
  
  return(p1) 
}

plotgg.site.diversity <- function(sitediv, p = c(0.025, 0.975),
                                  alpha = 0.2, theme = theme_blackbox, confints = TRUE){
  #sitediv <- sitediv.gen
  #p <- c(0.025,0.975)
  #alpha <- 0.2
  #theme <- theme_blackbox
  
  confint <- sitediv$mean + matrix(qnorm(p = c(0.025,0.975)),ncol = 2, nrow = nrow(sitediv),byrow = TRUE) * sitediv$sd
  colnames(confint) <- c("lower","upper")
  sitediv <- cbind(sitediv, confint)
  
  p1 <- ggplot(sitediv,aes(x = nsites, y = mean, col = group, group = group, fill = group)) +
    geom_line()
  
  if(confints){
    p1 <- p1 + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = alpha)
  }
  p1 <- p1 + theme
  #p1
  
  return(p1)
}


