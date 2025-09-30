#' Title
#'
#' @param ranks A data.table of variable used for ranking (decreasing= TRUE by default. See decreasing argument).
#' @param labels A data.table of logical labels (or that can be coerced to logical) for which ROC AUC will be computed .
#' @param compute.NES Should the Normalized Enrichment Score be computed (Niter= 1000)? Default= TRUE.
#' @param N.iter Number of iterations used to compute the NES. Default= 100L
#' @param decreasing Should the sort order for ranks be increasing or decreasing? Default= TRUE.
#' @param ranks.name The name to give to the ranks variables. Default= "rank".
#' @param labels.name The name to give to the labels variables. Default= "label".
#' @param cleanup.cache Should the temporary files be overwritten?
#'
#' @return
#' @export
#'
#' @examples
vl_rocAUC_grid <- function(ranks,
                           labels,
                           compute.NES= TRUE,
                           N.iter= 100L,
                           decreasing= TRUE,
                           ranks.name= "rank",
                           labels.name= "label",
                           cleanup.cache= FALSE)
{
  # Checks ----
  if(!is.data.table(ranks))
    ranks <- as.data.table(ranks)
  stopifnot(all(sapply(ranks, is.numeric)))
  if(!is.data.table(labels))
    labels <- as.data.table(labels)
  stopifnot(all(sapply(labels, function(x) is.numeric(x) | is.logical(x))))
  stopifnot(nrow(ranks)==nrow(labels))

  # Use a temp directory for caching signal files ----
  cache.file <- vl_cache_file(list(ranks, labels, compute.NES, N.iter, decreasing))

  # Main function ----
  if(cleanup.cache | !file.exists(cache.file)){

    # Loop over columns ----
    res <- lapply(
      seq_along(ranks),
      function(i) {
        message(paste0(i, "/", ncol(ranks), " rank"))
        lapply(
          seq_along(labels),
          function(j) {
            if(j %% 10 == 0)
              message(paste0(j, "/", ncol(labels), " label"))
            vl_rocAUC(
              ranks = ranks[[i]],
              labels = labels[[j]],
              compute.NES = compute.NES,
              N.iter = N.iter,
              decreasing = decreasing
            )
          }
        )
      }
    )

    # rbind result ----
    res <- lapply(res, function(x) {names(x) <- colnames(labels); rbindlist(x, idcol= labels.name)})
    names(res) <- colnames(ranks)
    res <- rbindlist(res, idcol= ranks.name)

    # Save ----
    saveRDS(res, cache.file)

  } else {

    # Import cached file ----
    res <- readRDS(cache.file)
  }

  # Format ----
  setattr(res, "class", c("vl_auc", "data.table", "data.frame"))

  # Return ----
  return(res)
}
