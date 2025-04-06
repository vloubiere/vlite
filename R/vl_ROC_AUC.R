#' ROC AUC
#'
#' @param predicted Predicted values from the model (ranging from 0 to 1).
#' @param label A vector of logical labels (or that can be coerced to logical).
#' @param plot Should the ROC be plotted? Default= FALSE.
#' @param xlab xlab. Default= "FALSE positive rate".
#' @param ylab ylab. Default= "TRUE positive rate".
#' @param type The type of plot, when add= FALSE. Default= "l".
#' @param add When plot is set to TRUE, should only the line be added to an existing plot?
#' @param ... Extra arguments to be passed to plot (when add= FALSE) or lines (when add= TRUE).
#'
#' @return ROC AUC
#' @export
vl_ROC_AUC <- function(predicted,
                       label,
                       plot= FALSE,
                       xlab= "False Positive Rate",
                       ylab= "True Positive Rate",
                       type= "l",
                       add= FALSE, ...)
{
  if(!is.logical(label))
    label <- as.logical(label)
  if(sum(label)==0)
    warning(paste0(length(label), "/", length(label), " labels are set to FALSE"))

  # Make data table ----
  dat <- data.table(label= label,
                    predicted= predicted)

  # Order ----
  setorderv(dat, "predicted", -1)

  # Compute AUC ----
  dat[, FPR:= cumsum(!label)/sum(!label)]
  dat[, TPR:= cumsum(label)/sum(label)]
  dat[, ROC_AUC:= c(0, sapply(seq(.N)[-1], function(i) (FPR[i] - FPR[i-1]) * (TPR[i] + TPR[i-1]) / 2))]

  # Plot ----
  if(plot)
  {
    if(add)
    {
      lines(dat$FPR,
            dat$TPR,
            ...)
    }else
    {
      lines(dat$FPR,
            dat$TPR,
            xlab= xlab,
            ylab= ylab,
            type= type,
            ...)
    }
  }

  # AUC function
  # AUC <- pROC::roc(AUC$label, AUC$score)$auc

  return(round(sum(dat$ROC_AUC), 2))
}
