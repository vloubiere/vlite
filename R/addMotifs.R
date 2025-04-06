#' Add motifs to enrichment plots
#'
#' @param DT DT object output from vl_motif_enrich or vl_motif_cl_enrich
#' @param cex.width expansion factor for motif widths
#' @param cex.height expansion factor for motif heights
#'
#' @examples
#' For vl_enr object
#' pl <- plot(vl_enr)
#' vl_add_motif(pl)
#'
#' For vl_enr_cl object
#' pl <- plot(vl_enr)
#' vl_add_motif(pl$DT)
#'
#' @export
addMotifs <- function(DT,
                      cex.width= 1,
                      cex.height= 1,
                      lwd= 0.1)
{
  # Extract PWMs
  DT <- unique(DT[, .(variable, name, y)])
  mats <- vl_Dmel_motifs_DB_full[DT, pwms_perc, on= "motif_ID==variable"]
  mats <- lapply(mats, TFBSTools::as.matrix)

  # Compute plotting coordinates
  ax.width <- max(strwidth(DT$name, cex= par("cex.axis")))+diff(grconvertX(c(0, par("mgp")[2]+0.5), "lines", "user"))
  coor <- addSeqLogo(pwm = mats,
                     x = par("usr")[1]-ax.width,
                     y = DT$y,
                     cex.width = cex.width,
                     cex.height = cex.height,
                     min_content= 0.05)

  # Plot
  coor[, segments(xleft,
                  ybottom,
                  xright,
                  ybottom,
                  xpd= T,
                  col= "grey60",
                  lwd= lwd*par("lwd"))]
  invisible(coor)
}
