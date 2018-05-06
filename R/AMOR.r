# Functions for general use in the metagenomics project.
# Most functions are for handling abundance tables.



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