#' Plots Rsq or PCC coeff
#'
#' @param pcc PCC value. If numeric, will be rounded, while characters will be printed as is
#' @param pos Position. Default= "topleft" (see x argument in ?legend for more detail).
#' @param digits Rounding digits. Default= 2.
#' @param bty Box around the legend? Default= FALSE.
#' @param ...
#'
#' @return Add PCC as a legend on aplot
#' @export
#'
#' @examples
addPcc <- function(pcc,
                   pos= "topleft",
                   bty= "n",
                   digits= 2,
                   ...)
{
  if(is.numeric(pcc))
    pcc <- round(pcc, digits) else
      stop("pcc should be numeric")

  legend(pos,
         legend= bquote(italic(r) == .(pcc)),
         bty= bty,
         ...)
}
