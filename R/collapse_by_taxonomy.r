collapse_by_taxonomy <- function(...) UseMethod("collapse_by_taxonomy")

collapse_by_taxonomy.default <- function(Tab,Tax,level=4,FUN=sum,sepchar=";"){
  # Match taxonomy and table rows
  #row.names(Tax) <- as.character(Tax$ID)
  
  if(class(Tab) != "matrix")
    stop("ERROR: Tab must be a matrix object",call.=TRUE)
  if(nrow(Tab) != nrow(Tax))
    stop("ERROR: Number of rows in Tab and Tax do not match",call.=TRUE)
  if(any(row.names(Tab) != row.names(Tax)))
    stop("ERROR: Row names for Tab and Tax do not match",call.=TRUE)
  
  tax <- get_tax_level(Tax,level=level,sepchar=sepchar)
  Tab.collapsed <- collapse_matrix(x=Tab,groups=tax,dim=1,FUN=FUN)
  return(Tab.collapsed)
}

collapse_by_taxonomy.Dataset <- function(Dat,level=4,FUN=sum,sepchar=";"){
  res <- collapse_by_taxonomy.default(Tab = Dat$Tab,
                                                Tax = Dat$Tax,
                                                level = level,
                                                FUN = FUN,
                                                sepchar = sepchar)
  
  if(length(Dat$Map) > 0){
    res <- create_dataset(Tab=res,Map=Dat$Map)
  }
  
  return(res)
}
