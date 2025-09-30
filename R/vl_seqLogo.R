#' Plot PWM seqlogo
#'
#' Transforms a PWM into an ICM before plotting the logo,
#'
#' @param pwm A PWM.
#' @param name Name of the motif to plot as a title.
#'
#' @examples
#'
#' @export
vl_seqLogo <- function(PWM,
                       name= NA)
{
  # Checks ----
  if(class(PWM)[1]=="PWMatrix") {
    name <- TFBSTools::name(PWM)
    PWM <- TFBSTools::as.matrix(PWM)
  }
  if(!is.matrix(PWM))
    stop("PWM should be a numeric matrix corresponding to a PWM.")

  # PWM to ICM ----
  ICM <- pwmToICM(PWM)
  ICM <- as.data.table(ICM, keep.rownames = "base")
  .c <- melt(ICM, "base", variable.name = "xleft", value.name = "height")
  .c[, xleft:= as.numeric(xleft)]
  setorderv(.c, "height", 1)
  .c[, ytop:= cumsum(height), xleft]

  # Plot ----
  plot(NA,
       xlim= c(0.5, ncol(ICM)+.5),
       ylim= c(0, max(.c$ytop)),
       xlab= "Position",
       ylab= "Bits",
       main= name)
  .c[, {
    plotDNAletter(base[1],
                  xleft[1],
                  ytop[1],
                  1,
                  height[1])
  }, .(base, xleft, height, ytop)]
  title(main= name)
}
