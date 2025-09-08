#' Title
#'
#' @param rankings
#' @param labels
#'
#' @return
#' @export
#'
#' @examples
vl_rankAUC <- function(rankings,
                       labels,
                       rank.name= "rank",
                       label.name= "motif",
                       cleanup.cache= FALSE)
{
  if(!is.data.table(labels))
    labels <- as.data.table(labels)
  if(!is.data.table(rankings))
    rankings <- as.data.table(rankings)
  stopifnot(all(sapply(rankings, is.numeric)))
  stopifnot(all(sapply(labels, function(x) is.numeric(x) | is.logical(x))))
  stopifnot(nrow(labels)==nrow(rankings))

  # Use a temp directory for caching signal files ----
  cache.file <- vl_cache_file(list(rankings, labels))

  # Main function ----
  res <- if(cleanup.cache | !file.exists(cache.file)){

    # Melt tables ----
    rankings <- data.table::melt(
      rankings,
      measure.vars = colnames(rankings),
      variable.name = "rank",
      value.name = "value.1"
    )
    labels <- data.table::melt(
      labels,
      measure.vars = colnames(labels),
      variable.name = "label",
      value.name = "value.2"
    )

    # Compute AUC ----
    rankings[, {
      labels[, {
        as.list(
          vl_rocAUC(
            ranking_score = value.1,
            labels = value.2
          )
        )
      }, label]
    }, rank]
  } else
    readRDS(cache.file)

  # Format ----
  setnames(res, "label", label.name)
  setnames(res, "rank", rank.name)
  setattr(res, "class", c("vl_auc", "data.table", "data.frame"))

  # Return ----
  return(res)
}
