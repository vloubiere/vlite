#' See most recent error file in a directory
#'
#' @param metadata Path to an excel metadata sheet
#' @param first.col.name If metadata in in .xlsx format, this character string will be used to detect the starting line. Default= "user".
#'
#' @export
importMetadataSheet <- function(metadata,
                                first.col.name= "user")
{
  meta <- if(grepl(".txt$", metadata))
  {
    return(fread(metadata))
  }else if(grepl(".rds$", metadata))
  {
    return(readRDS(metadata))
  }else if(grepl(".xlsx$", metadata))
  {
    sheet <- readxl::read_xlsx(metadata)
    # Slip lines before first column names
    sel.rows <- cumsum(c(sheet[[1]]==first.col.name, NA))>0
    if(sum(sel.rows, na.rm= TRUE)==0)
      stop("The 'user' column, which should be the first of the excel sheet, could not be found!")
    sheet <- sheet[(sel.rows),]
    meta <- as.data.table(sheet[-1,])
    names(meta) <- unlist(sheet[1,])
    return(meta)
  }else
    stop("metadata file extension should be one of '.txt', '.rds' or '.xlsx'!")
}
