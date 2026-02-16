#' Upload a File to Dropbox
#'
#' Uses the Dropbox Uploader script to upload a file or folder to a specified Dropbox path.
#'
#' @param local.path A character string specifying the local file or folder path to upload.
#' @param remote.path A character string specifying the destination path in Dropbox.
#'
#' @return
#' No return value. Executes the Dropbox Uploader script to upload the file or folder.
#'
#' @examples
#' \dontrun{
#' dropboxUpload("local/file.txt", "/remote/folder/file.txt")  # Upload a file
#' dropboxUpload("local/folder", "/remote/folder")  # Upload a folder
#' }
#'
#' @export
dropboxUpload <- function(local.path,
                          remote.path= "./")
{
  # cmd <- paste("sh /groups/stark/vloubiere/apps/dropbox/dropbox_uploader.sh upload", local.path, remote.path)
  cmd <- paste("bash /usr/local/bin/dropbox_uploader.sh upload", local.path, remote.path)
  system(cmd)
}

