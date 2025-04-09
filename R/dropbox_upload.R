#' Upload a File to Dropbox
#'
#' Uses the Dropbox Uploader script to upload a file or folder to a specified Dropbox path.
#'
#' @param local_path A character string specifying the local file or folder path to upload.
#' @param remote_path A character string specifying the destination path in Dropbox.
#'
#' @return
#' No return value. Executes the Dropbox Uploader script to upload the file or folder.
#'
#' @examples
#' \dontrun{
#' vl_dropbox_upload("local/file.txt", "/remote/folder/file.txt")  # Upload a file
#' vl_dropbox_upload("local/folder", "/remote/folder")  # Upload a folder
#' }
#'
#' @export
vl_dropbox_upload <- function(local_path,
                              remote_path)
{
  cmd <- paste("sh  /groups/stark/vloubiere/apps/dropbox/dropbox_uploader.sh upload", local_path, remote_path)
  system(cmd)
}

