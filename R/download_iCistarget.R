#' Download iCisTarget report
#'
#' Download an iCisTarget report, uncompresses it and stores it in a subfolder
#'
#' @param url URL where the reports is to be found
#' @param dir Directory where the report should be saved.
#' @param folder_name Name of the unziped subfolder. By default, the downloaded subfolder will be named using the iCisTarget job name.
#'
#' @export
download_iCistarget <- function(url,
                                dir,
                                folder_name)
{
  cmd <- paste0("cd ", normalizePath(dir), ";wget ", url, "; unzip archive.zip; rm archive.zip")
  system(cmd)
  if(missing(folder_name))
    folder_name <- unlist(fread(paste0(normalizePath(dir), "/icistarget/statistics.tbl"))[1, 1])
  final <- paste0(normalizePath(dir), "/", folder_name)
  cmd <- paste("mv", paste0(normalizePath(dir), "/icistarget"), final)
  system(cmd)
}
