# Input/Output functions of the AMOR package are here
read.pool <- function(file){
    # Function that takes a text file pooling
    # map and reads it into R.
    fc <- file(file)
    Pool <- strsplit(readLines(fc), "\t")
    close(fc)
    return(Pool)
}

read.am <- function(file,format="am",taxonomy=FALSE){
    # Wrapper for read table with the default options
    if(format == "qiime"){
      tab <- read.table(file,skip=1,sep="\t",comment.char="",header=T,row.names=1)
    }else if(format == "am"){
      tab <- read.table(file,sep="\t",comment.char="",header=T,row.names=1)
    }else{
      stop("Unknown format\n",call = TRUE)
    }
    
    if(!(taxonomy == FALSE)){
      tax <- data.frame(ID = row.names(tab),Taxonomy = tab[,taxonomy],row.names=row.names(tab))
      tab[,taxonomy] <- NULL
      
      return(create_dataset(Tab = tab, Tax = tax))
    }
    
    return(create_dataset(Tab = tab))
}

write.qiime <- function(...) UseMethod("write.qiime")

write.qiime.default <- function(Tab,file){
    # Function that takes an OTU table, and writes a file in
    # QIIME format with it.
    first <- "# QIIME v1.3.0 OTU table"
    header <- colnames(Tab)
    header <- c("#OTU ID", header)
    header <- paste(header,collapse="\t")
    fileConn <- file(file)
    writeLines(c(first,header),fileConn,sep="\n")
    close(fileConn)
    write.table(Tab,file=file,col.names=F,row.names=T,sep="\t",quote=F,append=T)
}

write.qiime.Dataset <- function(Dat,file){
  write.qiime.default(Tab = Dat$Tab, file = file)
}