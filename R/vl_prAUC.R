#' PR AUC
#'
#' Computes Precision-Recall AUC and optionally plots the PR curve.
#'
#' @param rank A numeric variable used for ranking (higher = more positive when decreasing=TRUE).
#' @param label A vector of logical labels (or that can be coerced to logical).
#' @param compute.NES Should a Normalized Enrichment Score be computed by label permutation? Default= TRUE.
#' @param N.iter Number of iterations used to compute the NES. Default= 100L.
#' @param decreasing Should the sort order for rank be decreasing? If FALSE, ranks are inverted for PRROC. Default= TRUE.
#' @param plot Should the PR curve be plotted? Default= FALSE.
#' @param main Main title. Default= "PR".
#' @param xlab xlab. Default= "Recall (TPR)".
#' @param ylab ylab. Default= "Precision (PPV)".
#' @param xlim Plot xlim. Default= c(0, 1).
#' @param ylim Plot ylim. Default= c(0, 1).
#' @param type The type of plot. Default= "l".
#' @param add Should the PR curve be added to an existing plot?
#' @param ... Extra arguments passed to plot (when plot=TRUE and add=FALSE) or lines (when add=TRUE).
#'
#' @return PR AUC, NES and/or plots the PR curve.
#' 
#' @examples
#' # Labels
#' label <- c(1,0,0,0,0, 1,0,0,0,0)  # 2 positives, 8 negatives (imbalanced)
#' 
#' # Case A: good ROC AUC AND good PR (positives ranked at the very top)
#' rank_good <- c(0.99, 0.20,0.19,0.18,0.17, 0.98, 0.16,0.15,0.14,0.13)
#' 
#' # Case B: good ROC AUC but bad PR (one negative is ranked above both positives)
#' # Positives are still above *most* negatives, so ROC AUC stays high, but early precision is poor, so PR suffers.
#' rank_bad <- c(0.90, 0.99,0.20,0.19,0.18, 0.89, 0.17,0.16,0.15,0.14)
#' 
#' # Compare roc and pr AUCs
#' vl_par(mfrow= c(2,2))
#' vl_rocAUC(rank= rank_good, label= label, plot= T)
#' vl_prAUC(rank= rank_good, label= label, plot= T)
#' vl_rocAUC(rank= rank_bad, label= label, plot= T)
#' vl_prAUC(rank= rank_bad, label= label, plot= T)
#' 
#' @export
vl_prAUC <- function(rank,
                     label,
                     compute.NES= TRUE,
                     N.iter= 100L,
                     decreasing= TRUE,
                     plot= FALSE,
                     main= "PR",
                     xlab= "Recall (TPR)",
                     ylab= "Precision (PPV)",
                     type= "l",
                     add= FALSE,
                     ...)
{
  # Checks ----
  if(!is.logical(label))
    label <- as.logical(label)
  stopifnot(is.logical(label))
  stopifnot(length(label) == length(rank))
  
  # PRROC expects higher score = more positive
  score <- if(decreasing) rank else -rank
  
  # Compute PR AUC ----
  PR <- PRROC::pr.curve(
    scores.class0 = score[label],   # positives (TRUE)
    scores.class1 = score[!label],  # negatives (FALSE)
    curve = TRUE
  )
  auc <- as.numeric(PR$auc.integral)
  
  # Compute NES ---
  if(compute.NES) {
    set.seed(123)
    rdm <- sapply(seq(N.iter), function(i) {
      perm <- sample(label)
      PRp <- PRROC::pr.curve(
        scores.class0 = score[perm],
        scores.class1 = score[!perm],
        curve = FALSE
      )
      as.numeric(PRp$auc.integral)
    })
    NES <- (auc - mean(rdm)) / sd(rdm)
  }
  
  # Plot ----
  if(add)
  {
    lines(PR$curve[,1],  # recall
          PR$curve[,2],  # precision
          col= "red",
          ...)
  } else if(plot)
  {
    plot(PR$curve[,1],   # recall
         PR$curve[,2],   # precision
         xlab= xlab,
         ylab= ylab,
         type= type,
         main= main,
         ...)
  }
  
  # Return ----
  var <- if(compute.NES)
    data.table(auc= auc, NES= NES) else
      data.table(auc= auc)
  return(var)
}