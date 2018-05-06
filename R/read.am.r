#' Read abundance matrix
#' 
#' Reads an abundance matrix into a Dataset object
#' 
#' @param file File path of the abundance matrix.
#' @param format Format of the abundance table file.
#' Either "am" or "qiime", see description for details.
#' @param taxonomy Indicates whether a taxonomy assignment
#' column is present in the abundance matrix. If FALSE,
#' no taxonomy is expected, if a string or integer are passed,
#' then the column with the corresponding name or number is
#' expected to have taxonomy information.
#' @param simplify Logical, indicating whether the abundance
#' table should be coerced to a numerical matrix
#' 
#' @details This function is a wrapper for read.table(), and
#' reads tab-delimited text files that contain abundance matrices.
#' It supports two formats: 1) am, this format contains a first row
#' of sample IDs and and a first colum of taxon IDs; and 2) qiime,
#' this format is the original QIIME format, where the first line
#' contains a version identification string, the second row sample
#' IDs and the first column (after second row) taxon IDs. The first
#' two lines start witht he '#' character.
#' 
#' @return If simplify = FALSE, it returns a Dataset object, where the value
#' of the Tax element is determined by the taxonomy argiment.
#' 
#' If simplify = TRUE, it returns a numeric abundance matrix. Any taxonomic
#' information is lost in this case.
#' 
#' @author Sur Herrera Paredes
#' 
#' @seealso \code{\link{create_dataset}}
#' 
#' @export
read.am <- function(file,format = "am", taxonomy = FALSE, simplify = FALSE){
  # Wrapper for read table with the default options
  # Read table in
  if(format == "qiime"){
    tab <- read.table(file,skip=1,sep="\t",comment.char="",header=T,row.names=1)
  }else if(format == "am"){
    tab <- read.table(file,sep="\t",comment.char="",header=T,row.names=1)
  }else{
    stop("Unknown format\n",call = TRUE)
  }
  
  # Separate taxonomy
  if(!(taxonomy == FALSE)){
    tax <- data.frame(ID = row.names(tab),
                      Taxonomy = tab[,taxonomy],
                      row.names=row.names(tab))
    tab[,taxonomy] <- NULL
    Dat <- create_dataset(Tab = tab, Tax = tax)
  }else{
    Dat <- create_dataset(Tab = tab)
  }
  
  # Extract abundance matrix if simplify is passed
  if(simplify){
    Dat <- Dat$Tab
  }
  
  return(Dat)
}