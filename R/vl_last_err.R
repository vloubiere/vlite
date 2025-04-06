#' See most recent error file in a directory
#'
#' @param dir URL where the reports is to be found
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
