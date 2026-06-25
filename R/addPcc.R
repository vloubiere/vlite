#' Plots Rsq or PCC coeff
#' 
#' Wrapper around vl_legend to add pcc to a plot.
#' 
#' @param x Numeric x-coordinate for legend position or a keyword. Default= par("usr")[2].
#' @param y Numeric y-coordinate for legend position. Default= par("usr")[4].
#' @param pcc PCC value. If numeric, will be rounded, while characters will be printed as is
#' @param digits Rounding digits. Default= 2.
#' @param legend Character vector of legend labels.
#' @param fill Vector of colors for legend key boxes.
#' @param x.adj If specified, x plotting position is shifted by x.adj*line.width. Default= 0.
#' @param y.adj If specified, y plotting positions is shifted by y.adj*line.height. Default= 0.
#' @param x.log Should the x value be logged? Default= par("xlog").
#' @param y.log Should the x value be logged? Default= par("ylog").
#' @param bty Character specifying box type around legend. Default is "n" (no box).
#' @param xpd Logical controlling clipping. Default is TRUE (allow legend outside plot region).
#' @param cex Numeric scaling factor for legend text. Default is 7/12 (≈0.58).
#' @param border Color for borders of legend boxes. Default is NA (no
#'   borders).
#' @param ... Additional arguments passed to legend().
#'
#' @return Add PCC as a legend on aplot
#' @export
#'
#' @examples
addPcc <- function(
    x= "topleft",
    y,
    pcc,
    digits= 2,
    fill,
    x.adj = 0,
    y.adj = 0,
    x.log = par("xlog"),
    y.log = par("ylog"),
    bty = "n",
    xpd = T,
    cex = 7/12,
    border = NA,
    ...)
{
  if(is.numeric(pcc))
    pcc <- round(pcc, digits) else
      stop("pcc should be numeric")
  
  vl_legend(
    x= x,
    y= y,
    legend= bquote(italic(r) == .(pcc)),
    fill= fill,
    x.adj = x.adj,
    y.adj = y.adj,
    x.log = x.log,
    y.log = x.log,
    bty = bty,
    xpd = xpd,
    cex = cex,
    border = border,
    ...
  )
}
