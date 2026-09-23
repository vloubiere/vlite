#' ROC AUC
#'
#' A wrapper around pROC::roc and pROC::auc that computes the ROC AUC, NES and/or plots the ROC curve.
#' 
#' @param rank A numeric variable used for ranking (decreasing= TRUE by default. See decreasing argument).
#' @param label A vector of logical labels (or that can be coerced to logical).
#' @param compute.NES Should the Normalized Enrichment Score be computed (Niter= 1000)? Default= TRUE.
#' @param N.iter Number of iterations used to compute the NES. Default= 100L
#' @param decreasing Should the sort order for rank be decreasing? Default= TRUE.
#' @param Ntop If specified, only the top N ranks will be used.
#' @param plot Should the ROC be plotted? Default= FALSE.
#' @param main Main title Default= "ROC".
#' @param xlab xlab. Default= "FPR".
#' @param ylab ylab. Default= "TPR".
#' @param type The type of plot. Default= "l".
#' @param add Should the ROC line be added to an existing plot?
#' @param ... Extra arguments to be passed to plot (when plot= TRUE) or lines (when plot= FALSE and add= TRUE).
#' 
#' @examples
#' # Example code
#' vl_rocAUC(rank= 10:1, label= rep(c(T, F), each= 5), plot= T)
#' 
#' @return ROC AUC, NES and/or the ROC curve.
#' 
#' @export
vl_rocAUC <- function(
    rank,
    label,
    compute.NES= TRUE,
    N.iter= 100L,
    seed= 123,
    decreasing= TRUE,
    Ntop,
    plot= FALSE,
    main= "ROC",
    xlab= "FPR",
    ylab= "TPR",
    type= "l",
    add= FALSE,
    col= "black",
    ...
)
{
  # Checks ----
  if(!is.logical(label))
    label <- as.logical(label)
  stopifnot(is.logical(label))
  stopifnot(any(label), any(!label))
  stopifnot(length(label)==length(rank))
  
  # Compute ROC auc ----
  ROC <- pROC::roc(
    response= label,
    predictor= rank,
    direction= ifelse(decreasing, "<", ">"),
    levels= c(F, T),
    na.rm= TRUE
  )
  auc <- as.numeric(pROC::auc(ROC))
  
  # Compute NES ---
  if(compute.NES) {
    set.seed(seed)
    rdm <- sapply(seq(N.iter), function(i) {
      rdm.ROC <- pROC::roc(
        response= sample(label),
        predictor= rank,
        direction= ifelse(decreasing, "<", ">"),
        levels= c(F, T)
      )
      as.numeric(pROC::auc(rdm.ROC))
    })
    NES <- (auc - mean(rdm)) / sd(rdm)
  } else
    NA <- as.numeric(NA)
  
  # Plot ----
  x.FPR <- 1-ROC$specificities # FPR= 1-specificity
  y.TPR <- ROC$sensitivities # TPR= sensitivity
  if(add)
  {
    lines(
      x.FPR,
      y.TPR,
      col= col,
      ...
    )
  } else if(plot)
  {
    plot(
      x.FPR,
      y.TPR,
      xlab= xlab,
      ylab= ylab,
      type= type,
      main= main,
      col= col,
      ...
    )
    abline(0, 1, lty = 2, col = "grey")
  }
  
  # Return ----
  var <- data.table(auc= auc, NES= NES, x.FPR= x.FPR, y.TPR= y.TPR)
  return(var)
}
