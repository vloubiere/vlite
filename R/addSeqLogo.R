#' plot seqlogo
#'
#' plot seqlogo from pwm matrix
#'
#' @param pwm List of pwm matrices
#' @param x x positions
#' @param y positions (centered)
#' @param pos eith 2 (left) of 4 (right)
#' @param cex.width width expansion factor applied before plotting motifs
#' @param cex.height height expansion factor applied before plotting motifs
#' @param add Should the pwm be plot on the top of opened device? Default= T
#'
#' @examples
#' pwm <- matrix(c(0.33,0.21,0.33,0.13,0.1,0.42,0.38,0.1,0.26,0.26,0.27,0.21,0,0.03,0.19,0.78,0.1,0.05,0.1,0.75,0.24,0.05,0.18,0.53,0.8,0.04,0,0.16,0.13,0.16,0.02,0.69,0.04,0.05,0.7,0.21,0.24,0.09,0.57,0.1,0.02,0.8,0.15,0.03,0.22,0.28,0.31,0.19,0.35,0.26,0.26,0.13,0.19,0.33,0.26,0.22), nrow= 4)
#' rownames(pwm) <- c("A", "C", "G", "T")
#' plot.new()
#' addSeqLogo(pwm)
#'
#' @export
addSeqLogo <- function(pwm,
                       x,
                       y,
                       pos= 2,
                       cex.width= 1,
                       cex.height= 1,
                       add= T,
                       min_content= 0)
{
  if(!pos %in% c(2,4))
    stop("Unsupported pos value. Use either 2 (left) or 4 (right)")
  if(is.matrix(pwm))
    pwm <- list(pwm)

  # Make object and index
  obj <- data.table(pwm, x, y, cex.width, cex.height)
  obj[, idx:= .I]

  # Width only depends on cex
  obj[, width:= strwidth("M", cex= cex.width), cex.width]

  # For each base, compute xleft, xright, ytop, ybottom
  pl <- obj[, {
    # Import PWM and melt
    .c <- as.data.table(pwm[[1]], keep.rownames= "base")
    .c <- melt(.c, id.vars = "base")
    # Compute motif content per column and normalize importance
    .c[, content:= sum(value*log2(value/c(0.25, 0.25, 0.25, 0.25)), na.rm= T), variable]
    .c[, norm:= value*(content/max(content))]
    # Remove flanks with little content (5% max)
    setorderv(.c, "variable")
    .c <- .c[min(which(norm>min_content)):max(which(norm>min_content))]
    # xleft depends on the pos (2 or 4)
    if(pos==4)
    {
      # Already correclty sorted earlier
      .c[, xleft:= x+((.GRP-1)*width), variable]
    }else if(pos==2)
    {
      setorderv(.c, "variable", -1)
      .c[, xleft:= x-(.GRP*width), variable]
    }
    # Rank from lowest to biggest importance -> inscreasing ytop pos
    setorderv(.c, "norm")
    .h <- strheight("M", cex= cex.height)
    .c[, c("height", "ytop"):= {
      heights <- norm*.h
      .(heights, (y-.h/2)+cumsum(heights))
    }, variable]
  }, .(idx, y, width)]

  # Plot
  pl[, plotDNAletter(base[1],
                     xleft[1],
                     ytop[1],
                     width[1],
                     height[1]), .(base, xleft, ytop, width, height)]

  # Return object containing limits of each motif
  invisible(pl[, .(xleft= min(xleft),
                   ybottom= min(ytop-height),
                   xright= max(xleft+width),
                   ytop= max(ytop)), .(idx)])
}
