#' Plot PFM
#'
#' @param pfm
#'
#' @return
#' @export
#'
#' @examples
vl_percLogo <- function(pfm) {
  # Compute position
  dat <- melt(as.data.table(pfm, keep.rownames = T), id.vars = "rn")
  setorderv(dat, "value")
  dat[, left:= .GRP-1, variable]
  dat[, ytop:= cumsum(value), variable]
  # Initiate plot
  plot(
    NA,
    type= "n",
    xlim= c(0, ncol(pfm)),
    ylim= c(0, max(dat$ytop, na.rm = T)),
    ylab= "Frequency"
  )
  # Plot
  dat[, {
    vlite::plotDNAletter(letter = rn, xleft = left, ytop = ytop, height = value, width = 1)
    print("")
  }, .(rn, left, ytop, value)]
}
