#' Plot contribution scores matrix
#'
#' @param contrib.object A contribution object, as outputed by ?importContrib()
#' @param enr And enrichment object, as outputed by ?contribEnrich() on the sequences of the contrib.object.
#' @param seqlvl A seqlvl present in in enr$sig.mot.coor$seqlvls.
#' @param min.count The minimum number of significant instances across (see enr$sig.inst). Default= 3L
#' @param best.by The group.by column for which only the instance with the best padjust should be returned.
#' @param sel If specified, only the motifs for which best.by in sel will be plotted.
#' @param xlab Default= "nt"
#' @param ylab Default= "Contribution"
#' @param xlim Default= sequence length.
#' @param ylim Default= range(contrib).
#'
#' @export
contribPlotLogo <- function(contrib.object,
                                 enr,
                                 seqlvl,
                                 min.count= 3L,
                                 best.by= "motif",
                                 sel,
                                 xlim,
                                 ylim,
                                 xlab= "nt",
                                 ylab= "Contribution")
{
  # Retrieve enriched motifs ----
  all.mots <- data.table::copy(enr)
  setnames(all.mots,
           old= best.by,
           new= "best.by")
  all.mots <- all.mots[, sig.mot.coor[[1]], .(best.by, sig.inst, log2OR, padj)]

  # Checks ----
  if(is.numeric(seqlvl) || length(seqlvl)>1)
    stop("seqlvl should be a character or a factor vector of length 1.")
  if(length(levels(all.mots$seqlvls)) != nrow(contrib.object))
    stop("The number of seqlvls in enr$sig.mot.coor should match the number of rows in contrib.object.")
  if(!seqlvl %in% levels(all.mots$seqlvls))
    stop("The seqlvl should be one of levels(enr$sig.mot.coor$seqlvls)")

  # Select sequence of interest ----
  mots <- all.mots[seqlvls==seqlvl & sig.inst>min.count, .SD[which.min(padj)], best.by]
  if(missing(sel))
    sel <- unique(mots$best.by)
  contrib <- contrib.object[which(levels(mots$seqlvls)==seqlvl)]
  contrib <- contrib[, .(base= unlist(tstrsplit(toupper(seq[[1]]), "")),
                         score= unlist(score))]

  # Plotting vars ----
  if(missing(xlim))
    xlim <- c(0, nrow(contrib))
  if(missing(ylim))
    ylim <- range(contrib$score)

  # Remove features outside xlim ----
  contrib[, xleft:= .I-1]
  contrib <- contrib[between(xleft, xlim[1], xlim[2])]
  mots <- mots[(start-1)<xlim[2] & end>xlim[1]]
  mots[(start-1)<xlim[1], start:= xlim[1]+1]
  mots[end>xlim[2], end:= xlim[2]]

  # Plotting ----
  plot(NA,
       xlim= xlim,
       ylim= ylim,
       xlab= xlab,
       ylab= ylab,
       frame= FALSE)
  contrib[, {
    plotDNAletter(letter = base,
                  xleft = xleft,
                  ytop= score,
                  width = 1,
                  height = score)
  }, (contrib)]

  # Add motif boxes ----
  mots[best.by %in% sel, {
    if(.N)
    {
      rect(xleft = start-1,
           ybottom = ylim[1],
           xright = end,
           ytop = ylim[2])
      text((start-1+end)/2,
           ylim[2],
           best.by,
           pos= 3,
           xpd= T)
    }
  }]

  # Return motifs
  return(mots)
}


