site_diversity <- function(Dat,factor,group,divfun = total_richness,nperm=20){
  #nperm <- 3
  #factor <- "Genotype"
  #group <- "Soil"
  #Dat <- Dat.rar
  #divfun <- total_richness
  
  #factor <- "accession"
  #group <- "Soil"
  #divfun <- total_richness
  #divfun <- divfun
  
  divfun <- match.fun(divfun)
  Dat.temp <- remove_samples(Dat, row.names(Dat$Map)[ Dat$Map[,factor] != group ])
  #Dat.temp$Tab <- Dat.temp$Tab > 0
  Perms <- matrix(ncol = nrow(Dat.temp$Map), nrow = nperm)
  for (i in 1:nperm){
    #i <- 1
    permutation <- sample(row.names(Dat.temp$Map))
    for (j in 1:length(permutation)){
      #j <- 2
      mat <- matrix(Dat.temp$Tab[,permutation[1:j]],ncol = j)
      S <- divfun(mat)
      Perms[i,j] <- S
    }
    
  }
  #Perms
  Res <- apply(Perms,2,function(vec) data.frame(mean = mean(vec), sd = sd(vec)))
  Res <- do.call(rbind,Res)
  Res$nsites <- 1:j
  Res$group <- group
  class(Res) <- c("site.diversity","data.frame")
  return(Res)
}

compare_site_diversity <- function(...) UseMethod("compare_site_diversity")
compare_site_diversity.default <- function(Tab, Map, factor, divfun = total_richness, nperm = 20){
  #Tab <- Dat.rar$Tab
  #Map <- Dat.rar$Map
  #factor <- "Genotype"
  #divfun <- total_richness
  #nperm <- 20
  
  groups <- levels(Map[,factor])
  Dat <- create_dataset(Tab = Tab, Map = Map)
  Res <- NULL
  for(group in groups){
    res <- site_diversity(Dat = Dat, factor = factor, group = group, divfun = divfun, nperm = nperm)
    Res <- rbind(Res,res)
  }
  
  class(Res) <- c("site.diversity","data.frame")
  return(Res)
}

compare_site_diversity.Dataset <- function(Dat, factor, divfun = total_richness, nperm = 20){  
  Res <- compare_site_diversity(Tab = Dat$Tab, Map = Dat$Map, factor = factor, divfun = divfun, nperm = nperm)
  Res$group <- factor(Res$group, levels = levels(Dat$Map[,factor]))
  return(Res)
}
