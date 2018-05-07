#' Get tax level
#' 
#' Processes a vector of txonomy strings and
#' returns a vector with the specified level
#' 
#' It splits each taxonomy string according to separator and
#' returns the specified level. If there are fewer levels
#' than the specified number, it returns the value 'unclassified'
#'
#' @param Tax A data.frame with a column named 'Taxonomy' which is
#' a character vector with taxonomy strings.
#' @param level The level to return
#' @param sepchar The separator of each level in the
#' taxonomy string
#'
#' @return A vector of taxonomy strings truncated at the
#' specified level
#' 
#' @author Sur Herrera Paredes
#' @export
#'
#' @examples
#' tax <- data.frame(ID = letters[1:4],
#'                   Taxonomy = c('a;b;c;d','a;b','A;B;C','A;B;C;D'))
#' get_tax_level(tax, level = 2)
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