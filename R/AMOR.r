# Functions for general use in the metagenomics project.
# Most functions are for handling abundance tables.

#library(vegan)
#library(labdsv,lib.loc="/nas02/home/s/u/sur/lib/R/")
#library(labdsv)

beta_diversity <- function(X,method="Bray-Curtis",sample_dim=2){
  X <- as.matrix(X)
  
  if(method == "Bray-Curtis"){
    if(sample_dim == 1){
      X.t <- t(X)
    }else if(sample_dim == 2){
      X.t <- X
    }else{
      stop("Invalid sample_dim",call=TRUE)
    }
    X.d <- apply(X,sample_dim,function(vec,X.t){  
      num <- colSums(abs(vec - X.t))
      den <- colSums(vec + X.t)
      d <- num/den  
      d
    },X.t=X.t)
    
    X.d <- as.dist(X.d)
  }else{
    stop("Unknown method",call=TRUE)
  }
  
  return(X.d)
}

communityComposition <- function(OTUtab,dist=NULL,method,col="grey",patch=19,cor=TRUE,main=""){
    # Method: pca or pco for the moment.

    #Define is matrix is the abundance table
    #or a distance table. Create input for method
    if(class(OTUtab) == "dist"){
        print("communityComposition: Input is distance object. Proceeding to evaluate directly.")
        Mat <- OTUtab
    }else{
        if(!is.null(dist)){
    	    Mat <- vegdist(t(OTUtab),method=dist)
        }else{
            Mat <- t(OTUtab)
        }
    }
    
    # Choose method and perform
    if(method == "pca"){
        # Make sure the input for PCA is a matrix
        if(class(Mat) == "dist"){
            Mat <- as.matrix(Mat)
        }
        Mat.p <- pca(Mat,cor=cor)
        plot.pca(Mat.p,pch=patch,col=col,title=main)
        return(Mat.p)
    }else if(method == "pco"){
        # Check the required parameters are passed. Either
        # I passed a dist object or I passed the dist argument
        # and transformed it.
        if(class(Mat) != "dist"){
            stop("communityComposition: You must specify a dist argument for PCoA",call.=TRUE)
        }
        Mat.p <- pco(Mat)
        plot.pco(Mat.p,pch=patch,col=col,title=main)
        return(Mat.p)
    }
}

extendRow <- function(row,map,n){
    #Function to extend an abundance table. Imported from MW_pipeline

    # Takes a named vector from and OTU table with the OTU
    # id in the end, together with a metadata table and converts
    # the vector into a dataframe where each element of the row
    # becomes a row itself with all the information about the
    # sample where it comes from
    
    # Get OTU name and get the vector of counts
    otu_name <- as.character(row[n])
    row <- row[-n]
    
    # Create unique names for the rows of the new data frame
    rowNames <- paste(otu_name,names(row),sep=".")
    
    # Create data frame
    Dat <- data.frame(otu_freq=as.numeric(row),otu_name=rep(otu_name,n-1))
    row.names(Dat) <- rowNames
    
    #Sanity check
    if(length(names(row)) != dim(map)[1]){
        stop("extendRow: The map file dimension does not coincide with the row dimension")
    }
    if(any(names(row) != row.names(map))){
        stop("extendRow: The sample names do not coincide between the map and the row")
    }
    
    #Add metadata columns
    Dat <- cbind(Dat,map)
    
    return(Dat)
}

extendTable <- function(Dat,map){
  # Function that takes two data frames (an OTU table and a metadata
  # file), and generates and extended table.
  if(any(names(Dat) == "otu")){
    stop("extendTable: Invalid sample name (otu)")  
  }
  if(class(Dat) != "data.frame"){
    stop("extendTable: OTU table must be a data frame")
  }
  
  # Add otu id's to the data frame
  Dat$otu <- row.names(Dat)
  
  # length of OTU vector
  n <- dim(Dat)[2]
  
  # Extend row by row
  Profile <- do.call(rbind,apply(Dat,1,extendRow,map=map,n=n))
  return(Profile)
}

findGoodSamples <- function(OTUtab,min_reads_sample){
    # Finds OTUs with reads above threshold
    index <- colSums(OTUtab) > min_reads_sample
    return(index)
}

makeMultipleRarefactionGroups <- function(id,depth,names,rarefaction_threshold=2000){
    # This is mainly a wrapper for when you want to pool multiple groups
    # of samples, each one independently
    Meta <- data.frame(id,depth,names)
    Meta$id <- as.character(Meta$id)
    Res <- NULL
    for(name in levels(Meta$names)){
        Dat <- Meta[ Meta$names == name, ]
        Res <- c(Res,makeRarefactionGroups(id=as.character(Dat$id),name=name,depth=Dat$depth,rarefaction_threshold=rarefaction_threshold))
    }
    return(Res)
}

makeRarefactionGroups <- function(name,id,depth,rarefaction_threshold=2000){
    # This function takes a set of samples and the sampling depth
    # on each, and a descriptor, and generates and object that specifies
    # how to pool the samples to get as many pooled samples above
    # the specified threshold.
    Group <- data.frame(id,depth)
    Group$id <- as.character(Group$id)
    Group <- Group[order(Group$depth,decreasing=TRUE),]
    #Flag for the first loop
    continue <- TRUE
    #Indexes move in opposite direction until they collapse
    pos1 <- 0
    pos2 <- length(Group$depth)
    #Variable that tracks the depth of the current group
    curr_depth <- 0
    Pools <- NULL
    while(continue){
        pos1 <- pos1 + 1
        curr_depth <- Group$depth[pos1]
        curr_group <- Group$id[pos1]
        flag <- TRUE
        while(flag){
            #First check the indexes to figure out where we are
            if(pos2 < pos1){
                #If it is smaller, it means the sample was already added
                flag <- FALSE
                continue <- FALSE
            }else if(pos2 == pos1){
                # If they are the same it means the current samples haven't
                # been added. Merge with prev group
                if(!is.null(Pools)){
                    Pools[[length(Pools)]] <- c(Pools[[length(Pools)]],curr_group)
                }
                flag <- FALSE
                continue <- FALSE            
            } else if(curr_depth + Group$depth[pos2] >= rarefaction_threshold){
                # Finish group
                curr_group <- c(curr_group,Group$id[pos2])
                Pools[[length(Pools) + 1]] <- c(name,curr_group)
                #restart
                pos2 <- pos2 - 1
                prev_group <- curr_group
                curr_group <- NULL
                curr_depth <- 0
                flag <- FALSE
            }else{
                #updating depth and samples in group
                curr_depth <- curr_depth + Group$depth[pos2]
                curr_group <- c(curr_group,Group$id[pos2])
                #updating pos2 index
                pos2 <- pos2 - 1
            }
        }
    }
    #if(is.null(Pools)){
    #    print(id)
    #}
    return(Pools)
}

normalizeSample <- function(sample){
    100*sample/sum(sample)
}

summarizeOTUdistribution <- function(OTUtab){
  # Eventually should show how OTUs distribute
  # in samples.
  print("Not implemented")
}


##### EXTERNAL FUNCTIONS

get_tax_level <- function(Tax,level=4,sepchar=";"){
  tax <- strsplit(as.character(Tax$Taxonomy),split=sepchar)
  tax <- sapply(tax,function(vector,level){
    if(length(vector) >= level){
      paste(vector[1:level],collapse=sepchar)
    }else{
      "unclassified"
    }},level=level)
  
  tax
}