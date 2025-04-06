#' PR AUC
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
#' @return PR AUC
#' @export
vl_PR_AUC <- function(predicted,
                      label,
                      plot= FALSE,
                      xlab= "True Positive Rate",
                      ylab= "Positive Predictive Value",
                      type= "l",
                      add= FALSE, ...)
{
  if(!is.logical(label))
    label <- as.logical(label)
  if(sum(label)==0)
    warning(paste0(length(label), "/", length(label), " labels are set to FALSE"))

  # Make data table ----
  dat <- data.table(label= label,
                    pred= predicted)

  # Order ----
  setorderv(dat, "pred", -1)

  # Compute AUC ----
  dat[, TPR:= cumsum(label)/sum(label)]
  dat[, PPV:= cumsum(label)/seq(.N)]
  dat[, PR_AUC:= c(0, sapply(seq(.N)[-1], function(i) (TPR[i] - TPR[i-1]) * (PPV[i] + PPV[i-1]) / 2))]

  # Plot ----
  if(plot)
  {
    if(add)
    {
      lines(dat$TPR,
            dat$PPV,
            ...)
    }else
    {
      plot(dat$TPR,
           dat$PPV,
           xlab= xlab,
           ylab= ylab,
           type= type,
           ...)
    }
  }

  # AUC function
  # AUC <- pROC::roc(AUC$label, AUC$score)$auc

  return(round(sum(dat$PR_AUC), 2))
}
