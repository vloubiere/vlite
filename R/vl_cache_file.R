#' Title
#'
#' @param input.list list of elements used by the function which will be used to generate a reproducible path for the cache file.
#' @param suffix A prefix to be added at the beginning of the cache file name. Default= "tmp_".
#' @param extension Extension for the cache file. Default= ".rds".
#' @param tmp.dir Default= tempdir().
#'
#' @return
#' @export
#'
#' @examples
vl_cache_file <- function(input.list,
                          prefix= "tmp_",
                          extension= ".rds",
                          tmp.dir= tempdir()) {
  # Check
  if(is.vector(input.list))
    input.list <- as.list(input.list)
  stopifnot(is.list(input.list))
  if(!is.null(prefix))
    prefix <- paste0(prefix, "_")

  # Use a temp directory for caching signal files ----
  signal.key <- digest::digest(input.list)
  cache.path <- file.path(tmp.dir, paste0(prefix, signal.key, extension))

  # Return cache path
  return(cache.path)
}
