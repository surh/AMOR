phylogram <- function(...) UseMethod("phylogram")

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

phylogram.Dataset <- function(Dat, facet = NULL,colname="Sample",
                              variable.name="Taxon",value.name="Abundance",
                              scales="free_x",space = "free_x", nrow.legend=20,
                              ntaxa = NULL,other_name = "other"){
  res <- phylogram.default(Tab = Dat$Tab, Map = Dat$Map, facet = facet, colname = colname,
                    variable.name = variable.name, value.name = value.name, scales = scales,
                    space = space, nrow.legend = nrow.legend, ntaxa = ntaxa,
                    other_name = other_name)
  
  return(res)
}
