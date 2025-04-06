#' TRUE positive rate curve
#'
#' @param predicted Predicted values from the model (ranging from 0 to 1).
#' @param label A vector of logical labels (or that can be coerced to logical).
#' @param Nleft Number of enhancers left before cutoff. Default= 100.
#' @param plot Should the TPR be plotted?
#' @param xlim x limits for plotting. Default= NULL.
#' @param ylim y limits for plotting. Default= NULL.
#' @param xlab x label. Default= "Prediction score".
#' @param ylab y label. Default= "Positive pred. value (percentage)".
#' @param lty.1 Line type before Nleft. Default= 1.
#' @param lty.2 Line type after Nleft. Default= 3.
#' @param col Color.
#' @param cex.max cex for the max TPR text label.
#' @param pos.max pos for the TPR text label Default= 3.
#' @param offset.max Offset for plotting the TPR text label. Default= .5.
#' @param add When plot is set to TRUE, should only the lines/label be added to an existing plot?
#' @param ... Extra arguments to be passed to lines
#'
#' @return TPR plot
#' @export
vl_TPR <- function(predicted,
                   label,
                   Nleft= 100,
                   plot= FALSE,
                   xlim= NULL,
                   ylim= NULL,
                   xlab= "Prediction score",
                   ylab= "TRUE positive rate (%)",
                   lty.1= 1,
                   lty.2= 3,
                   col= "black",
                   cex.max= .7,
                   pos.max= 3,
                   offset.max= .5,
                   add= FALSE,
                   ...)
{
  if(!is.logical(label))
    label <- as.logical(label)
  if(sum(label)==0)
    warning(paste0(length(label), "/", length(label), " labels are set to FALSE"))

  # Create a data table with observed and predicted values
  dat <- data.table(label = label,
                    pred = predicted)

  # Calculate the cumulative percentage of actually positive values
  setorderv(dat, cols = "pred") # Sort predicted values in ascending order
  dat[, FN:= cumsum(label)] # FALSE negative
  setorderv(dat, cols = "pred", order = -1) # Sort predicted values in descending order
  dat[, TP:= cumsum(label)] # TRUE positive
  dat[, TPR:= TP/(TP+FN)*100]
  setorderv(dat,
            cols = "pred")

  # Fit a smooth spline and find peaks
  x <- dat[1:(.N-Nleft), pred]
  y <- dat[1:(.N-Nleft), TPR]
  spline_fit <- smooth.spline(x, y)
  spline_derivative <- predict(spline_fit, deriv = 1)

  # Find peaks
  peaks <- which(diff(sign(spline_derivative$y)) == -2) + 1
  peak_x_values <- spline_derivative$x[peaks]
  peak_y_values <- predict(spline_fit, x = peak_x_values)$y
  peak_x_values <- c(peak_x_values, last(x)) # value at .N-Nleft
  peak_y_values <- c(peak_y_values, last(y)) # value at .N-Nleft

  # Select ideal cutoff
  sel <- min(which(peak_y_values>max(0.95*peak_y_values)))
  x_cutoff <- peak_x_values[sel]
  y_cutoff <- peak_y_values[sel]

  # Compute limits
  if(is.null(xlim))
    xlim <- range(dat$pred)
  if(is.null(ylim))
    ylim <- c(0, min(c(100, max(dat$TPR)*1.1)))

  # Plot TPR and cutoffs
  if(plot)
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

    # Plot TPR lines
    dat[,{
      # TPR
      lines(pred[1:(.N-Nleft)],
            TPR[1:(.N-Nleft)],
            lty= lty.1,
            col= col,
            ...)
      lines(pred[(.N-Nleft):.N],
            TPR[(.N-Nleft):.N],
            lty = lty.2,
            col= col,
            ...)
    }]
  }

  # Return cutoffs
  return(list(min_TPR= dat$TPR[1],
              predict_cutoff= x_cutoff,
              TPR_at_cutoff= y_cutoff))
}
