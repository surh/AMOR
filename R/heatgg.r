heatgg <- function(...) UseMethod("heatgg")

heatgg.default <- function(Tab, Map, order.samples.by = NULL, facet = NULL, sample.id.var = "SAMPLEID",
    col.name = "Taxon", value.name = "Abundance", trans = "log",
    guide = "colourbar", gradientn.colours = c("white","#67001F"),
                           facet.scales = "free",cluster = FALSE,
                           distfun = dist){
  # thanks to
  # https://github.com/chr1swallace/random-functions/blob/master/R/ggplot-heatmap.R
  # http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2
  
  # Check for errors
  if(!all(colnames(Tab) == row.names(Map))){
    stop("ERROR: Samples names in Tab do not math sample names in Map.", call. = TRUE)
  }
  if(sample.id.var %in% colnames(Map)){
    stop("ERROR: sample.id.var already exists in Map.",call. = TRUE)
  }
  if(!is.null(order.samples.by)){
    if(!order.samples.by %in% colnames(Map)){
      stop("ERROR: order.samples.by does not exist in Map", call. = TRUE)
    }
  }
  if(!is.null(facet) & cluster){
    stop("ERROR: faceting and clustering are not supported at the same time.",call. = TRUE)
  }
  
  # Combine abundance and metadata
  id.vars <- c(colnames(Map),sample.id.var)
  Dat <- cbind(t(Tab), Map)
  Dat[,sample.id.var] <- row.names(Dat)
  
  # Cluster if needed
  if(cluster){
    # Cluster
    taxon.hc <-hclust(distfun(Tab))
    sample.hc <- hclust(distfun(t(Tab)))
    taxon.dd <- dendro_data(as.dendrogram(taxon.hc))
    sample.dd <- dendro_data(as.dendrogram(sample.hc))
    
    # Plot dendograms
    taxon.dd.plot <- ggplot(segment(taxon.dd)) +
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      labs(x = NULL, y = NULL) +
      theme_dendro() +
      theme(plot.background = element_blank(),
            plot.margin = unit(c(0,0,0,0),"mm"))
    #taxon.dd.plot
      
    sample.dd.plot <- ggplot(segment(sample.dd)) +
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
      scale_x_continuous(breaks = 1:ncol(Tab),labels=sample.hc$labels[sample.hc$order], expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme(panel.background = element_blank(),
            plot.background = element_blank(),
            axis.text.y = element_text(angle = 0),
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.margin = unit(c(0,0,0,0),"mm")) +
      coord_flip()
      #scale_y_reverse()
    #sample.dd.plot
    
    # Get orders
    sample.ord <- colnames(Tab)[ match(sample.dd$labels$label,colnames(Tab)) ]
    taxon.ord <- colnames(t(Tab))[ match(taxon.dd$labels$label,row.names(Tab)) ]
    
    # Order taxon and make sur eno reordering happens
    #Tab <- Tab[ taxon.ord, sample.ord]
    #Map <- Map[ sample.ord, ]
    order.samples.by <- NULL
  }
  
  # Sort
  if(!is.null(order.samples.by)){
    sample_order <- Dat[ order(Dat[,order.samples.by]), sample.id.var ]
    Dat[,sample.id.var] <- factor(Dat[,sample.id.var], levels = sample_order)
  }
  
  # Prepare for plotting
  Dat2 <- melt(Dat,id.vars=id.vars,variable.name = col.name, value.name = value.name)
  
  # Reorder for clustering
  if(cluster){
    Dat2$Taxon <- factor(Dat2$Taxon, levels = taxon.ord)
    Dat2[,sample.id.var] <- factor(Dat2[,sample.id.var], levels = sample.ord)
  }
  
  # Plot heatmap tiles
  p1 <- ggplot(Dat2,aes_string(x=col.name,y=sample.id.var))
  if(!is.null(facet)){
    p1 <- p1 + facet_grid(facets = facet, scales = facet.scales)
  }
  p1 <- p1 + geom_tile(aes_string(fill=value.name)) +
    scale_fill_gradientn(colours= gradientn.colours,guide = guide,trans = trans,na.value = NA) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme(panel.background = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(face = "bold",colour="black",size=7),
          axis.text.x = element_blank(),
          legend.position = "right",
          legend.key.size = unit(0.02,units = "npc"))

  if(cluster){
    # Combine dendrograms with tile plot
    p1 <- list(p1 = p1, taxon.dd.plot = taxon.dd.plot, sample.dd.plot = sample.dd.plot)
    class(p1) <- "heatggclus"
    #print(p1)
  }
  
  return(p1)
}

heatgg.Dataset <- function(Dat, order.samples.by = NULL, facet = NULL, sample.id.var = "SAMPLEID",
                           col.name = "Taxon", value.name = "Abundance", trans = "log",
                           guide = "colourbar", gradientn.colours = c("white","#67001F"),
                           facet.scales = "free", cluster = FALSE, distfun = dist){
  p1 <- heatgg.default(Tab = Dat$Tab, Map = Dat$Map,order.samples.by = order.samples.by,
                 facet = facet,sample.id.var = sample.id.var, col.name = col.name,
                 value.name = value.name, trans = trans, guide = guide,
                 gradientn.colours = gradientn.colours, facet.scales = facet.scales, cluster = cluster,
                 distfun = distfun)
  return(p1)
}
  
#### UTILITIES ###
print.heatggclus <- function(x,row.width = 0.2, col.width = 0.2){
  grid.newpage()
  top.layout <- grid.layout(nrow = 2, ncol = 2,
                            widths = unit(c(1-row.width,row.width), "null"),
                            heights = unit(c(col.width,1-row.width), "null"))
  pushViewport(viewport(layout=top.layout))
  if(col.width>0)
    print(x$taxon.dd.plot, vp=viewport(layout.pos.col=1, layout.pos.row=1))
  if(row.width>0)
    print(x$sample.dd.plot, vp=viewport(layout.pos.col=2, layout.pos.row=2), width = row.width)
  ## print centre without legend
  print(x$p1 + labs(x = NULL,y = NULL) +
          theme(axis.ticks = element_blank(),
                axis.line = element_blank(),
                axis.text = element_blank(),
                axis.title.x=element_blank(),axis.title.y=element_blank(),
                legend.position="none",
                panel.border=element_blank(),panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),plot.background=element_blank(),
                #plot.margin = unit(c(0,0,0,0),"mm")),
                plot.margin = unit(rep(0,4),"lines")),
        vp=viewport(layout.pos.col=1, layout.pos.row=2))
  
  ## add legend
  legend <- g_legend(x$p1)
  pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
  grid.draw(legend)
  upViewport(0)
  
  res <- list(p1 = x$p1, taxon.dd.plot = x$taxon.dd.plot, sample.dd.plot = x$sample.dd.plot)
  class(res) <- "heatggclus"
  return(res)
}

g_legend <- function(a.gplot){
  # http://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

