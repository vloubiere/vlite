#' plot seqlogo
#'
#' plot seqlogo from pwm matrix
#'
#' @param pwm A percentage PWM where colsums= 1 (or PWMatrix).
#' @param name Name of the motif (will be used as title).
#'
#' @examples
#' pwm <- matrix(
#' c(0.33,0.21,0.33,0.13,0.1,0.42,0.38,0.1,0.26,0.26,0.27,0.21,0,0.03,
#' 0.19,0.78,0.1,0.05,0.1,0.75,0.24,0.05,0.18,0.53,0.8,0.04,0,0.16,
#' 0.13,0.16,0.02,0.69,0.04,0.05,0.7,0.21,0.24,0.09,0.57,0.1,0.02,0.8,
#' 0.15,0.03,0.22,0.28,0.31,0.19,0.35,0.26,0.26,0.13,0.19,0.33,0.26,0.22
#' ), nrow= 4)
#' rownames(pwm) <- c("A", "C", "G", "T")
#' plot.new()
#' addSeqLogo(pwm)
#'
#' @export
vl_seqLogo <- function(pwm,
                       name= NA)
{
  # Checks ----
  if(class(pwm)=="PWMatrix") {
    name <- TFBSTools::name(pwm)
    pwm <- TFBSTools::as.matrix(pwm)
  }
  if(!is.matrix(pwm))
    stop("pwm should be a PWMatrix or a matrix object.")
  if(is.null(colnames(pwm)))
    colnames(pwm) <- seq(ncol(pwm))

  # Format data.table object ----
  .c <- as.data.table(pwm, keep.rownames= "base")
  .c <- melt(.c, id.vars = "base")

  # Compute Information content and height
  .c[, IC := 2 + sum(ifelse(value > 0, value * log2(value), 0)), by = variable]
  .c[, height := value * IC]

  # Compute plotting parameters ----
  setorderv(.c, "height") # Highest scores on top
  .c[, xleft:= as.numeric(variable)-.5]
  .c[, ytop:= cumsum(height), variable]
  .c[, width:= 1]

  # Plot ----
  plot(NA,
       xlim= c(0.5, ncol(pwm)+.5),
       ylim= c(0, max(.c[, sum(height), variable]$V1)),
       xlab= "Position",
       ylab= "Weighted contribution")
  .c[, {
    plotDNAletter(base[1],
                  xleft[1],
                  ytop[1],
                  width[1],
                  height[1])
  }, (.c)]
  title(main= name)
}
