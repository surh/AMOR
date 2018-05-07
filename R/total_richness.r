#' Calculate richness of an abundance table.
#' 
#' Given a numeric matri representign an abundance
#' table, it calculates the total richness (ie
#' species number) in it.
#'
#' @param Tab A numeric matrix with columns as samples and rows as taxa.
#' 
#' @return A number representing the total richness.
#' @export
#' @author Sur Herrera Paredes
#' @seealso \link{site_diversity}, \link{compare_site_diversity}
total_richness <- function(Tab){
  sum(rowSums(Tab) > 0)
}