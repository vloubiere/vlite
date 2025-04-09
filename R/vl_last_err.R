#' View the Most Recent Error File
#'
#' Displays the most recent `.err` file in a specified directory.
#'
#' @param dir A character string specifying the folder path where `.err` files are stored.
#'
#' @return
#' No return value. Opens the most recent `.err` file for viewing and prints its modification time.
#'
#' @examples
#' \dontrun{
#' vl_last_err("/path/to/error/files")  # View the most recent .err file in the directory
#' }
#'
#' @export
vl_last_err <- function(dir)
{
  dat <- data.table(file= list.files(dir, ".err$", full.names = T))
  dat[, mtime:= file.info(file)$mtime]
  last <- dat[which.max(mtime)]
  print(paste("Printed on:", last$mtime))
  file.show(last$file)
}
