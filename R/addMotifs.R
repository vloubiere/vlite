#' Add motifs to enrichment plots
#'
#' @param plot.DT Output of the plot.vl_enr or plot.vl_enr_cl methods containing plotting coordinates.
#' @param pwms A list of 'PWMatrix' whose names correspond to the 'motif' entries in plot.DT.
#' @param cex.width expansion factor for motif width.
#' @param cex.height expansion factor for motif height.
#'
#' @examples
#' For vl_enr object
#' pl <- plot(vl_enr)
#' vl_add_motif(pl)
#'
#' For vl_enr_cl object
#' pl <- plot(vl_enr)
#' vl_add_motif(pl$plot.DT)
#'
#' @export
addMotifs <- function(plot.DT,
                      pwms,
                      cex.width= 1,
                      cex.height= 1,
                      lwd= 0.1)
{
  # Extract pwms ----
  idx <- match(plot.DT$motif, sapply(pwms, TFBSTools::name))

  # Check ----
  if(anyNA(idx))
    stop("Some motifs in plot.DT$motif could not be found in sapply(pwms, TFBSTools::name)")

  # To matrix ----
  mats <- lapply(pwms[idx], TFBSTools::as.matrix)

  # Compute plotting coordinates ----
  lab.width <- strwidth(plot.DT$name, cex= par("cex.axis"))
  space.plot.lab <- diff(grconvertX(c(0, par("mgp")[2]+0.5), "lines", "user"))
  lab.left <- max(lab.width)+space.plot.lab
  mot.right <- par("usr")[1]-lab.left

  # Plot motifs ----
  coor <- addSeqLogo(pwm = mats,
                     x = mot.right,
                     y = plot.DT$y,
                     cex.width = cex.width,
                     cex.height = cex.height,
                     min_content= 0.05)

  # Add lines under motifs ----
  coor[, segments(xleft,
                  ybottom,
                  xright,
                  ybottom,
                  xpd= T,
                  col= "grey60",
                  lwd= lwd*par("lwd"))]
}
