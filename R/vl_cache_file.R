#' Title
#'
#' @param input.list list of elements used by the function which will be used to generate a reproducible path for the cache file.
#'
#' @return
#' @export
#'
#' @examples
vl_cache_file <- function(input.list) {
  # Check
  if(is.vector(input.list))
    input.list <- as.list(input.list)
  stopifnot(is.list(input.list))

  # Use a temp directory for caching signal files ----
  cache_dir <- tempdir()
  signal.key <- digest::digest(input.list)
  cache.path <- file.path(cache_dir, paste0(signal.key, ".rds"))

  # Return cache path
  return(cache.path)
}
