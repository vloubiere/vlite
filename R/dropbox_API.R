#' Use dropbox uploader
#' @export
vl_dropbox_upload <- function(local_path,
                              remote_path)
{
  cmd <- paste("sh  /groups/stark/vloubiere/apps/dropbox/dropbox_uploader.sh upload", local_path, remote_path)
  system(cmd)
}

#' Use dropbox downloader
#' @export
vl_dropbox_download <- function(local_path,
                                remote_path)
{
  cmd <- paste("sh  /groups/stark/vloubiere/apps/dropbox/dropbox_uploader.sh download", remote_path, local_path)
  system(cmd)
}
