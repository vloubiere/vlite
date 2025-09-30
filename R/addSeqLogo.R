#' Add ICM logos to an existing plot
#'
#' plot seqlogo from pwm matrix
#'
#' @param pwm List of percentage PWM matrices which will be converted to ICM before plotting.
#' @param x x positions (see pos argument).
#' @param y y positions (centered).
#' @param pos Either 2 (left) of 4 (right)
#' @param cex.width width expansion factor.
#' @param cex.height height expansion factor.
#' @param min.content Flanks with a smaller summed content are not plotted.
#'
#' @examples
#' # Example: start from a PPM and create a PWM, then plot ICM
#' PPM <- matrix(
#'   c(
#'     0.33,0.21,0.33,0.13,0.1,0.42,0.38,0.1,0.26,0.26,0.27,0.21,0,0.03,
#'     0.19,0.78,0.1,0.05,0.1,0.75,0.24,0.05,0.18,0.53,0.8,0.04,0,0.16,
#'     0.13,0.16,0.02,0.69,0.04,0.05,0.7,0.21,0.24,0.09,0.57,0.1,0.02,0.8,
#'     0.15,0.03,0.22,0.28,0.31,0.19,0.35,0.26,0.26,0.13,0.19,0.33,0.26,0.22
#'   ), nrow= 4
#' )
#' PWM <- freqToPWM(PPM)
#'
#' # Plot
#' plot(0, 1, type= "n")
#' addSeqLogo(PWM, x = -1, y= 1, cex.height = 3, cex.width = 2, pos= 2, min_content = .5)
#'
#' @export
addSeqLogo <- function(pwm,
                       x,
                       y,
                       pos= 2,
                       cex.width= 1,
                       cex.height= 1,
                       min.content= 0.05)
{
  # Checks ----
  if(!pos %in% c(2,4))
    stop("Unsupported pos value. Use either 2 (left) or 4 (right)")
  if(is.matrix(pwm))
    pwm <- list(pwm)
  # Cast TFBS objects to matrix
  class <- unique(sapply(pwm, function(x) class(x)[1]))
  if(class %in% c("PWMatrix", "PFMatrix"))
    pwm <- lapply(pwm, TFBSTools::as.matrix)

  # Compute width and height ----
  .w <- strwidth("M", cex= cex.width)
  .h <- strheight("M", cex= cex.height)

  # Loop ----
  lapply(seq_along(pwm), function(i) {
    # Convert pwm to ICM
    ICM <- pwmToICM(pwm[[i]])
    # Mean content filtering
    sel <- colSums(ICM)>min.content
    ICM <- ICM[, min(which(sel)):max(which(sel))]
    # Cast
    ICM <- as.data.table(ICM, keep.rownames = "base")
    # Compute coordinates
    .c <- melt(ICM, "base", variable.name = "xleft", value.name = "height")
    .c[, xleft:= x[i]+as.numeric(xleft)*.w]
    setorderv(.c, "height", 1)
    .c[, height:= height/max(height)*.h]
    .c[, ytop:= y[i]+cumsum(height)-.h/2, xleft]
    # Revert
    if(pos==4)
      .c[, xleft:= xleft-(ncol(ICM)+2)*.w]
    # Plot
    .c[, {
      plotDNAletter(base[1],
                    xleft[1],
                    ytop[1],
                    .w,
                    height[1])
    }, .(base, xleft, height, ytop)]
  }
  )
}
