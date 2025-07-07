#' Scatterplot function
#'
#' A wrapper around ?plot with nicer default settings.
#'
#' @param x
#' @param y Default= NULL
#' @param pch Default= 19
#' @param cex Default= .8
#' @param col Default= adjustcolor("grey", .7)
#' @param xaxt Default= "s"
#'
#' @return
#' @export
#'
#' @examples
#' vl_par()
#' vl_plot(1:3)
vl_plot <- function(x,
                    y= NULL,
                    pch= 19,
                    cex= .8,
                    col= adjustcolor("grey", .7),
                    xaxt= "s",
                    frame= F,
                    ...)
{
  plot(x= x,
       y= y,
       pch= pch,
       cex= cex,
       col= col,
       xaxt= "n",
       frame= frame,
       ...)
  if(xaxt!="n")
    axis(1,
         padj= -1.25)
}
