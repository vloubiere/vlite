#' Positive Predicted Value
#'
#' Given a set of predicted values and a vector of logical labels, computes the Positive Predictive Values for
#' increasing predicted thresholds.
#'
#' @param predicted Predicted values from the model (ranging from 0 to 1).
#' @param label A vector of logical labels (or that can be coerced to logical).
#' @param Nleft Minimum number of remaining positive sequences (prediction>=cutoff) before computing max PPV.
#' Default= 100.
#' @param plot Should the PPV be plotted?
#' @param xlim x limits for plotting. Default= NULL.
#' @param ylim y limits for plotting. Default= NULL.
#' @param xlab x label. Default= "Prediction score".
#' @param ylab y label. Default= "Positive pred. value (percentage)".
#' @param lty.1 Line type before Nleft. Default= 1.
#' @param lty.2 Line type after Nleft. Default= 3.
#' @param col Color.
#' @param cex.max cex for the max PPV text label.
#' @param pos.max pos for the PPV text label Default= 3.
#' @param offset.max Offset for plotting the PPV text label. Default= .5.
#' @param add When plot is set to TRUE, should only the lines/label be added to an existing plot?
#' @param return.value.at If specified, return the PPV value at the last predicted value <= return.value.at.
#' @param ... Extra arguments to be passed to lines
#'
#' @return PPV at a given threshold and/or PPV curve with max PPV with a minimum number of observations.
#' 
#' @examples
#' # Create synthetic data
#' predicted <- c(1,1,0.8,0.4,.1,0)
#' label <- c(T,T,T,F,F,F)
#' 
#' # Plot PPV
#' vl_PPV(predicted, label, Nleft= 3, plot= T) # 100% TPR with at least 3 obs
#' vl_PPV(predicted, label, Nleft= 3, return.value.at= .4) # 75% TPR with pred >= .4 (and at least 3 obs)  
#' vl_PPV(predicted, label, Nleft= 4, plot= T) # 100% TPR with at least 4 obs
#' 
#' @export
vl_PPV <- function(predicted,
                   label,
                   Nleft= 100,
                   plot= FALSE,
                   xlim= NULL,
                   ylim= NULL,
                   xlab= "Prediction score",
                   ylab= "Positive pred. value (%)",
                   lty.1= 1,
                   lty.2= 3,
                   col= "black",
                   show.max= TRUE,
                   show.max.value= show.max,
                   cex.max= .7,
                   pos.max= 2,
                   offset.max= .5,
                   show.pred.cutoff= T,
                   add= FALSE,
                   return.value.at= NULL,
                   ...)
{
  # Checks ----
  if(!is.logical(label))
    label <- as.logical(label)
  if(!is.null(return.value.at) && !is.numeric(return.value.at))
    warning("return.value.at should be numeric")
  if(plot && !is.nul(return.value.at))
    message("When return.value.at is specified, plot is ignored.")
  
  # Create a data table with observed and predicted values
  dat <- data.table(label = label,
                    pred = predicted)
  if(anyNA(dat)) {
    warning("NA values will be dropped!")
    dat <- na.omit(dat)
  }
  
  # Extra checks ----
  if(sum(dat$label)==0)
    warning("All labels are FALSE")
  stopifnot(Nleft <= nrow(dat))
  
  # Sort the data table by the predicted values in descending order
  setorderv(dat, cols = "pred", order = -1)
  
  # Aggregate for each predicted value (handle rank ties)
  agg <- dat[, .(`TRUE.lab`= sum(label), `FALSE.lab`= sum(!label)), by= pred]
  
  # Calculate the cumulative percentage of actually positive values
  agg[, TP:= cumsum(TRUE.lab)] # TRUE positive
  agg[, FP:= cumsum(FALSE.lab)] # FALSE positive
  agg[, PPV:= TP/(TP+FP)*100]
  setorderv(agg, "pred")
  
  # Return value at certain cutoff
  if(!is.null(return.value.at)) {
    # If any predicted value passes the threshold
    PPV.cutoff <- if(any(agg$pred >= return.value.at)) {
      # Return 
      agg[pred >= return.value.at][1, PPV]
    } else {
      # No value
      as.numeric(NA)
    }
    return(PPV.cutoff)
  }
  
  # Find max value
  max.PPV <- which.max(agg[(TP+FP)>=Nleft, PPV])
  x_cutoff <- agg$pred[max.PPV]
  y_cutoff <- agg$PPV[max.PPV]
  
  # Compute limits
  if(is.null(xlim))
    xlim <- range(agg$pred)
  if(is.null(ylim))
    ylim <- c(0, min(c(100, max(agg$PPV)*1.1)))
  
  # Plot PPV and cutoffs
  if(plot | add)
  {
    # Initiate plot
    if(!add)
    {
      plot(NA,
           NA,
           type = "n",
           xlab = xlab,
           ylab = ylab,
           xlim= xlim,
           ylim= ylim,
           ...)
    }
    
    # Plot PPV lines
    agg[,{
      # PPV
      lines(pred[(TP+FP)>=Nleft],
            PPV[(TP+FP)>=Nleft],
            lty= lty.1,
            col= col,
            ...)
      lines(pred[(TP+FP)<=Nleft],
            PPV[(TP+FP)<=Nleft],
            lty = lty.2,
            col= col,
            ...)
      
      # Plot top PPV point
      if(show.max) {
        points(x_cutoff,
               y_cutoff,
               col = adjustcolor(col, .7),
               pch = 19)
      }
      
      # Plot top PPV value
      if(show.max.value) {
        text(x_cutoff,
             y_cutoff,
             paste0(round(y_cutoff, 1), "%"),
             pos= pos.max,
             offset= offset.max,
             col= col,
             cex= cex.max)
      }
      
      
      # Plot x prediction cutoff
      if(show.pred.cutoff) {
        segments(x_cutoff,
                 0,
                 x_cutoff,
                 y_cutoff,
                 lty= "33")
        text(x_cutoff,
             0,
             round(x_cutoff, 2),
             pos= 2,
             offset= offset.max,
             cex= cex.max)
      }
    }]
  }
  
  # Return cutoffs
  return(list(min_PPV= agg$PPV[1],
              predict_cutoff= x_cutoff,
              PPV_at_cutoff= y_cutoff))
}
