#' ROC AUC
#'
#' @param ranks A vector of numeric variables used for ranking (decreasing= TRUE by default. See decreasing argument).
#' @param labels A vector of logical labels to compute the ROC (or that can be coerced to logical).
#' @param compute.NES Should the Normalized Enrichment Score be computed (Niter= 1000)? Default= TRUE.
#' @param N.iter Number of iterations used to compute the NES. Default= 100L
#' @param decreasing Should the sort order for ranks be increasing or decreasing? Default= TRUE.
#' @param plot Should the ROC be plotted? Default= FALSE.
#' @param main Main title Default= "ROC".
#' @param xlab xlab. Default= "FPR".
#' @param ylab ylab. Default= "TPR".
#' @param type The type of plot. Default= "l".
#' @param add Should the ROC line be added to an existing plot?
#' @param ... Extra arguments to be passed to plot (when plot= TRUE) or lines (when plot= FALSE and add= TRUE).
#'
#' @return ROC AUC
#' @export
vl_rocAUC <- function(ranks,
                      labels,
                      compute.NES= TRUE,
                      N.iter= 100L,
                      decreasing= TRUE,
                      plot= FALSE,
                      main= "ROC",
                      xlab= "FPR",
                      ylab= "TPR",
                      type= "l",
                      add= FALSE,
                      ...)
{
  # Checks ----
  if(!is.logical(labels)) {
    labels <- as.logical(labels)
    warning("labels coerced to logical")
  }
  stopifnot(length(ranks) == length(labels))

  # Compute N pos and neg ----
  n_pos <- sum(labels)
  n_neg <- sum(!labels)
  stopifnot(n_pos > 0, n_neg > 0)

  # Order ----
  labels <-  labels[order(ranks, decreasing= decreasing, method= "radix")]

  # Compute TRUE and FALSE positive rate ----
  TPR <- c(0, cumsum(labels) / n_pos)
  FPR <- c(0, cumsum(!labels) / n_neg)

  # AUC via trapezoidal rule ----
  auc <- sum(diff(FPR) * (head(TPR, -1) + tail(TPR, -1)) / 2)

  # bg for NES ---
  set.seed(123)
  if(compute.NES) {
    rdm <- sapply(seq(N.iter), function(i) {
      perm <- sample(labels)
      pTPR <- c(0, cumsum(perm) / n_pos)
      pFPR <- c(0, cumsum(!perm) / n_neg)
      sum(diff(pFPR) * (head(pTPR, -1) + tail(pTPR, -1)) / 2)
    })
    NES <- (auc - mean(rdm)) / sd(rdm)
  }

  # Plot ----
  if(add)
  {
    lines(FPR,
          TPR,
          ...)
  } else if(plot)
  {
    plot(FPR,
         TPR,
         xlab= xlab,
         ylab= ylab,
         type= type,
         main= main,
         ...)
    abline(0, 1, lty = 2, col = "grey")
  }

  # Return ----
  var <- if(compute.NES)
    data.table(auc= auc, NES= NES) else
      data.table(auc= auc)
  return(var)
}
