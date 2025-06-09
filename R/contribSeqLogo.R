#' Plot contribution scores matrix
#'
#' @param contrib A contrib object imported using ?importContrib.
#' @param row.idx The line index within the contrib object that should be plotted. Default= 1.
#' @param start Integer specifying at which the start position of the sequence should be clipped.
#' @param end Integer specifying at which the end position of the sequence should be clipped.
#' @param higlight Index of the positions that should be higlighted.
#'
#' @export
contribSeqLogo <- function(contrib,
                           row.idx= 1,
                           start,
                           end,
                           highlight.pos)
{
  # Get sequence and contrib ----
  .c <- contrib[row.idx]
  base <- unlist(tstrsplit(toupper(.c$seq), ""))
  score <- unlist(.c$score)

  # Clip sequence if necessary ----
  if(missing(start))
    start <- 1
  if(missing(end))
    end <- length(base)
  base <- base[start:end]
  score <- score[start:end]
  highlight <- highlight.pos[between(highlight.pos, start, end)]

  # Initiate plot ----
  plot(NA,
       xlim= c(start-1, end),
       ylim= range(score),
       xlab= "Position",
       ylab= "Contribution score")
  segments(par("usr")[1], 0, par("usr")[2], 0, lty= 3)

  # Highlight rectangles ----
  breaks <- c(TRUE, diff(highlight) != 1, TRUE)
  starts <- highlight[breaks[-length(breaks)]]
  ends <- highlight[breaks[-1]]
  rect(xleft = starts,
       ybottom = par("usr")[3],
       xright = ends,
       ytop = par("usr")[4],
       col= "lightgrey",
       border= NA)

  # Contrib ----
  sapply(seq(base), function(i) {
    .s <- score[i]
    plotDNAletter(letter = base[i],
                  xleft = start+i-1,
                  ytop= ifelse(.s<0, 0, .s),
                  width = 1,
                  height = ifelse(.s<0, -.s, .s))
  })

  invisible()
}


