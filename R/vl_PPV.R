#' Positive Predicted Value curve
#'
#' @param predicted Predicted values from the model (ranging from 0 to 1).
#' @param label A vector of logical labels (or that can be coerced to logical).
#' @param Nleft Number of enhancers left before cutoff. Default= 100.
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
#' @param ... Extra arguments to be passed to lines
#'
#' @return PPV plot
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

  # Sort the data table by the predicted values in descending order
  setorderv(dat, cols = "pred", order = -1)

  # Calculate the cumulative percentage of actually positive values
  dat[, TP:= cumsum(label)] # TRUE positive
  dat[, FP:= cumsum(!label)] # FALSE positive
  dat[, PPV:= TP/(TP+FP)*100]
  setorderv(dat,
            cols = "pred")

  # Fit a smooth spline and find peaks
  x <- dat[1:(.N-Nleft), pred]
  y <- dat[1:(.N-Nleft), PPV]
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
    ylim <- c(0, min(c(100, max(dat$PPV)*1.1)))

  # Plot PPV and cutoffs
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

    # Plot PPV lines
    dat[,{
      # PPV
      lines(pred[1:(.N-Nleft)],
            PPV[1:(.N-Nleft)],
            lty= lty.1,
            col= col,
            ...)
      lines(pred[(.N-Nleft):.N],
            PPV[(.N-Nleft):.N],
            lty = lty.2,
            col= col,
            ...)

      # Plot the smooth spline
      # lines(spline_fit, col = "blue", lwd = 2)

      # Plot top PPV point
      points(x_cutoff,
             y_cutoff,
             col = adjustcolor(col, .7),
             pch = 19)
      text(x_cutoff,
           y_cutoff,
           paste0(round(y_cutoff, 1), "%"),
           pos= pos.max,
           offset= offset.max,
           col= col,
           cex= cex.max)

      # Add segments
      # segments(x_cutoff,
      #          0,
      #          x_cutoff,
      #          y_cutoff,
      #          lty= "33")
      # segments(0,
      #          y_cutoff,
      #          x_cutoff,
      #          y_cutoff,
      #          lty= "33")
      # text(x_cutoff,
      #      y_cutoff/2,
      #      round(x_cutoff, 2),
      #      pos= 4,
      #      offset= 1)
      # text(x_cutoff/2,
      #      y_cutoff,
      #      paste0(round(y_cutoff, 1), "%"),
      #      pos= 3,
      #      offset= 1)
    }]
  }

  # Return cutoffs
  return(list(min_PPV= dat$PPV[1],
              predict_cutoff= x_cutoff,
              PPV_at_cutoff= y_cutoff))
}
