#' Plots Rsq or PCC coeff
#'
#' @param rsq Rsquare value. If numeric, will be rounded, while characters will be printed as is
#' @param pos Position. Default= "topleft" (see x argument in ?legend for more detail).
#' @param digits Rounding digits. Default= 2.
#' @param adjusted Is the Rsq adjusted? Default= FALSE.
#' @param bty Box around the legend? Default= FALSE.
#' @param ...
#'
#' @return Add Rsquare as a legend on aplot
#' @export
#'
#' @examples
addRsq <- function(rsq,
                   pos= "topleft",
                   adjusted= F,
                   bty= "n",
                   digits= 2,
                   ...)
{
  if(is.numeric(rsq))
    rsq <- round(rsq, digits) else
      stop("rsq should be numeric")

  if(adjusted)
  {
    legend(pos,
           legend= bquote(Adj.~R^2 == .(rsq)),
           bty= bty,
           ...)
  }else
  {
    legend(pos,
           legend= bquote(R^2 == .(rsq)),
           bty= bty,
           ...)
  }
}
