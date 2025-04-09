#' Download a File or Folder from Dropbox
#'
#' Uses the Dropbox Uploader script to download a file or folder from a specified Dropbox path to a local path.
#'
#' @param local_path A character string specifying the local destination path for the downloaded file or folder.
#' @param remote_path A character string specifying the source path in Dropbox.
#'
#' @return
#' No return value. Executes the Dropbox Uploader script to download the file or folder.
#'
#' @examples
#' \dontrun{
#' vl_dropbox_download("local/file.txt", "/remote/folder/file.txt")  # Download a file
#' vl_dropbox_download("local/folder", "/remote/folder")  # Download a folder
#' }
#'
#' @export
vl_dropbox_download <- function(local_path,
                                remote_path)
{
  cmd <- paste("sh  /groups/stark/vloubiere/apps/dropbox/dropbox_uploader.sh download", remote_path, local_path)
  system(cmd)
}
