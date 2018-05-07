#' Write a QIIME abundance table file
#' 
#' Writes a file compatible with QIIME
#'
#' @param x Either an abundance matrix or a Dataset
#' @param file Path to the file to write
#'
#' @export
#' @author Sur Herrera Paredes
#'
#' @examples
#' data(Rhizo)
#' 
#' # The following are equivalent
#' write.qiime(Rhizo,'myfile.txt')
#' write.qiime(create_dataset(Rhizo),'myfile.txt')
write.qiime <- function(x, file) UseMethod("write.qiime")

#' @rdname write.qiime
#' @method write.qiime default
write.qiime.default <- function(x, file){
    # Function that takes an OTU table, and writes a file in
    # QIIME format with it.
    first <- "# QIIME v1.3.0 OTU table"
    header <- colnames(x)
    header <- c("#OTU ID", header)
    header <- paste(header, collapse="\t")
    fileConn <- file(file)
    writeLines(c(first, header), fileConn,sep="\n")
    close(fileConn)
    write.table(x, file=file, col.names=F,
                row.names=T, sep="\t",
                quote=F, append=T)
}

#' @rdname write.qiime
#' @method write.qiime Dataset
write.qiime.Dataset <- function(x, file){
  write.qiime.default(x = x$Tab, file = file)
}