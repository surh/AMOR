#' Collapse by taxonomy
#' 
#' Collapses taxa according to the taxonomy object or a given factor
#' 
#' @param Tab An abundance matrix with samples as columns and taxa as rows
#' @param Tax A Taxonomy data.frame as defined by \link{create_dataset}
#' @param Group If NULL, then the "Taxonomy" column of the Taxonomy data.frame will be used.
#' If not NULL, it can be either a column number or name (i.e. character string), indicating
#' which column of the Taxonomy data.frame to use as grouping factor. It can also be a factor
#' or vector that can be coerced into factor.
#' @param level If Group is NULL, the taxonomy level to use from the Taxonomy column.
#' @param FUN Which function to use to aggregate
#' @param sepchar Separator character for taxonomy levels.
#' 
#' @author Sur Herrera Paredes
#' 
#' @rdname collapse_by_taxonomy
#' @export collapse_by_taxonomy
collapse_by_taxonomy <- function(...) UseMethod("collapse_by_taxonomy")

#' @rdname collapse_by_taxonomy
#' @method collapse_by_taxonomy default
#' @S3method collapse_by_taxonomy default
collapse_by_taxonomy.default <- function(Tab,Tax,Group = NULL, level=4,FUN=sum,sepchar=";"){
  # Match taxonomy and table rows
  #row.names(Tax) <- as.character(Tax$ID)
  
  if(class(Tab) != "matrix")
    stop("ERROR: Tab must be a matrix object",call.=TRUE)
  if(nrow(Tab) != nrow(Tax))
    stop("ERROR: Number of rows in Tab and Tax do not match",call.=TRUE)
  if(any(row.names(Tab) != row.names(Tax)))
    stop("ERROR: Row names for Tab and Tax do not match",call.=TRUE)
  
  if(is.null(Group) && is.numeric(level)){
    tax <- get_tax_level(Tax,level=level,sepchar=sepchar)
  }else if(!is.null(Group)){
    if(length(Group) == 1){
      tax <- Tax[ ,Group ]
    }else{
      tax <- factor(Group)
    }
  }else{
    stop("ERROR: Either a non-null group must be provided, or a numeric level for taxonomyu collapsint",
         call. = TRUE)
  }
  
  
  Tab.collapsed <- collapse_matrix(x=Tab,groups=tax,dim=1,FUN=FUN)
  return(Tab.collapsed)
}

#' @rdname collapse_by_taxonomy
#' @method collapse_by_taxonomy Dataset
#' @S3method collapse_by_taxonomy Dataset
collapse_by_taxonomy.Dataset <- function(Dat,Group, level=4,FUN=sum,sepchar=";"){
  res <- collapse_by_taxonomy.default(Tab = Dat$Tab,
                                      Tax = Dat$Tax,
                                      Group = Group,
                                      level = level,
                                      FUN = FUN,
                                      sepchar = sepchar)
  
  if(length(Dat$Map) > 0){
    res <- create_dataset(Tab=res,Map=Dat$Map)
  }
  
  return(res)
}
