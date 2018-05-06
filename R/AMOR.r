# Functions for general use in the metagenomics project.
# Most functions are for handling abundance tables.

normalizeSample <- function(sample){
    100*sample/sum(sample)
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