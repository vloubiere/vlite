#' Plot pval
#'
#' Plots formated pval data as text
#'
#' @param x X position for plotting
#' @param y Y position for plotting
#' @param pval Pval to be plotted
#' @param stars Should stars be printed?
#' @param values Should p values be printed?
#' @param pos pos argument (label position). see ?text()
#' @param offset P values offest (default= 0)
#' @param cex cex expansion factor passed to text
#' @param ... Extra arguments passed to test function
#' @examples
#' pval <- c(1e-10, 1e-5, 1e-3, 1e-2, 1e-1, 1)
#' plot(NA, xlim= c(0,7), ylim= c(-1,1))
#' vl_plot_pval(1:6, 0.5, pval)
#' @export
addPval <- function(x,
                    y,
                    pval,
                    stars= T,
                    values= F,
                    pos= 3,
                    offset= ifelse(values, -.2, -.35),
                    cex= .6,
                    ...)
{
  # Compute Stars ----
  star <- if(stars)
    cut(pval,
        breaks = c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf),
        labels = c("****", "***", "**", "*", "N.S")) else
          rep("", length(pval))
  star <- as.character(star)

  # Format pval ----
  lab <- if(values)
    formatC(pval, digits = 1, format= "e") else
      rep("", length(pval))

  # Compute labels ----
  mapply(function(x, y, pos, offset, cex, p, l, s)
  {
    var <- if(values)
    {
      if(p<2.2e-308)
        bquote(italic(P) < "2.2e-308" * .(s)) else if(p>0.05)
          bquote(italic(P) == .(l)^"N.S") else
            bquote(italic(P) == .(l) * .(s))
    }else
    {
      if(p>0.05)
        bquote(.(l)^"N.S") else
          bquote(.(l) * .(s))
    }
    text(x,
         y,
         labels= var,
         offset= offset,
         pos= pos,
         cex= cex,
         xpd= NA,
         ...)
  }, x= x, y= y, pos= pos, offset= offset, cex= cex, p= pval, l= lab, s= star)
}
