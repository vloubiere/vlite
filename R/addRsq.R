#' Plots Rsq or PCC coeff
#'
#' @param x Rsquare value. If numeric, will be rounded, while characters will be printed as is
#' @param pos Position. Default= "topleft" (see x argument in ?legend for more detail).
#' @param digits Rounding digits. Default= 2.
#' @param adjusted Is the Rsq adjusted? Default= FALSE.
#' @param bty Box around the legend? Default= FALSE.
#' @param ...
#'
#' @return Add Rsquare x as a legend on aplot
#' @export
#'
#' @examples
addRsq <- function(x,
                   pos= "topleft",
                   type= "rsq",
                   adjusted= F,
                   bty= "n",
                   digits= 2,
                   ...)
{
  if(is.numeric(x))
    x <- round(x, digits) else
      stop("x should be numeric")

  if(adjusted)
  {
    legend(pos,
           legend= bquote(Adj.~R^2 == .(x)),
           bty= bty,
           ...)
  }else
  {
    legend(pos,
           legend= bquote(R^2 == .(x)),
           bty= bty,
           ...)
  }
}
