#' Download a File or Folder from Dropbox
#'
#' Uses the Dropbox Uploader script to download a file or folder from a specified Dropbox path to a local path.
#'
#' @param remote.path A character string specifying the source path in Dropbox.
#' @param local.path A character string specifying the local destination path for the downloaded file or folder.
#'
#' @return
#' No return value. Executes the Dropbox Uploader script to download the file or folder.
#'
#' @examples
#' \dontrun{
#' dropboxDownload("/remote/folder/file.txt", "local/file.txt")  # Download a file
#' dropboxDownload("/remote/folder", "local/folder")  # Download a folder
#' }
#'
#' @export
dropboxDownload <- function(remote.path,
                            local.path)
{
  if(!dir.exists(dirname(local.path)))
    warning("local directory does not exist!")
  # cmd <- paste("sh /groups/stark/vloubiere/apps/dropbox/dropbox_uploader.sh download", remote.path, local.path)
  cmd <- paste("bash /usr/local/bin/dropbox_uploader.sh download", remote.path, local.path)
  system(cmd)
}
