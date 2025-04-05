#' Barplot allowing to add points and sd
#'
#' @param height A vector of values describing the bars which make up the plot.
#' @param individual.var List of individual variables to be plotted as points.
#' @param show.sd Should sd from individual var be shown as whiskers around bars?
#' @param bar.top.labels Labels to be shown on the top or each bar.
#' @param bar.top.labels.cex Expansion factor for bar.top.labels.
#' @param xlim
#' @param ylim
#' @param sd.arrow.length Length of sd arrow length. defaults to 1/8 of distance between bar centers.
#' @param sd.arrow.lwd Line width of sd arrows Default= .5.
#' @param individual.var.pch pch value for individual variables
#' @param individual.var.col color value for individual variables
#' @param individual.var.jitter jitter value for individual variables
#' @param individual.var.cex cex factor for points
#' @param compute.bar.diff List of pairwise bar pairs to be compared (V2/V11 ratio)
#' @param compute.bar.diff.digits Number of compute.bar.diff.digits for the computed diff. Default= 2
#' @param compute.bar.diff.cex cex value for pairwise comparisons
#' @param xpd xpd value for individual.var, sd, pairwaise comparisons...
#' @param horiz Not supported atm
#' @param ... Extra arguments to be passed to barplot
#'
#' @return A simple barplot
#' @export
#'
#' @examples
#' vl_barplot(1:3, rep(.2, 3))
vl_barplot <- function(height,
                       individual.var= NULL,
                       show.sd= !is.null(individual.var),
                       bar.top.labels= NULL,
                       compute.bar.diff= NULL,
                       xlim= NULL,
                       ylim= NULL,
                       tilt.names= TRUE,
                       individual.var.pch= 16,
                       individual.var.col= adjustcolor("lightgrey", .7),
                       individual.var.jitter= .2,
                       individual.var.cex= .5,
                       sd.arrow.lwd= .5,
                       sd.arrow.length= ifelse(length(height)>1, diff(bar[c(1,2)])/8, .2),
                       bar.top.labels.cex= .7,
                       compute.bar.diff.digits= 2,
                       compute.bar.diff.cex= .8,
                       xpd= NA,
                       horiz= F,
                       xaxt= "s",
                       names.arg= names(height),
                       ...)
{
  # Checks
  if(is.table(height))
    height <- c(height)
  if(!is.vector(height))
    stop("Only height tables and vectors are supported for now")
  if(horiz)
    stop("horiz not supported yet ;)")

  # Compute sd if necessary ----
  sd <- if(show.sd && !is.null(individual.var))
    sapply(individual.var, sd, na.rm= T) else
      rep(0, length(height))

  # Create data table to keep track of all values ----
  dat <- data.table(height= height,
                    sd.min= height-sd,
                    sd.max= height+sd)
  dat[, min:= apply(.SD, 1, function(x) min(unlist(x), na.rm= T))]
  dat[, max:= apply(.SD, 1, function(x) max(unlist(x), na.rm= T))]

  # Compute ylim that takes into account outliers ----
  if(is.null(ylim))
  {
    ylim <- range(c(dat$min, dat$max), na.rm = T)
    if(ylim[1]>0)
      ylim[1] <- 0
    if(ylim[2]<0)
      ylim[2] <- 0
  }

  # Initiate barplot ----
  bar <- barplot(height,
                 xlim= xlim,
                 ylim= ylim,
                 xaxt= ifelse(tilt.names, "n", xaxt),
                 names.arg= names.arg, ...)
  if(tilt.names && !is.null(names.arg) && xaxt!="n")
    tiltAxis(bar, labels= names.arg)
  dat[, x:= bar]

  # Add sd arrows if specified ----
  if(show.sd)
  {
    segments(bar,
             height-sd,
             bar,
             height+sd,
             lwd= sd.arrow.lwd,
             xpd= xpd)
    segments(bar-sd.arrow.length,
             height+sd,
             bar+sd.arrow.length,
             height+sd,
             lwd= sd.arrow.lwd,
             xpd= xpd)
    segments(bar-sd.arrow.length,
             height-sd,
             bar+sd.arrow.length,
             height-sd,
             lwd= sd.arrow.lwd,
             xpd= xpd)
  }

  # Add individual measurements if specified ----
  if(!is.null(individual.var))
  {
    if(!is.list(individual.var))
      individual.var <- as.list(individual.var)
    if(length(individual.var)!=length(height))
      stop("individual.var should be a vector or a list of the same length as height")
    x <- rep(bar, lengths(individual.var))
    y <- unlist(individual.var)
    points(jitter(x, amount = individual.var.jitter),
           y,
           pch= individual.var.pch,
           col= unlist(individual.var.col),
           xpd= xpd,
           cex= individual.var.cex)
  }

  # Add mean diff if specified ----
  if(!is.null(compute.bar.diff))
  {
    # Pairs to compare
    comp <- do.call(rbind, compute.bar.diff)
    comp <- data.table(V1= comp[,1],
                       V2= comp[,2],
                       var1= dat$height[comp[,1]],
                       var2= dat$height[comp[,2]],
                       x0= bar[comp[,1]],
                       x1= bar[comp[,2]])
    comp[, x:= rowMeans(.SD), .SDcols= c("x0", "x1")]
    # Compute FC and max/ values (y position)
    comp[, var:= paste0("x", formatC(var2/var1, compute.bar.diff.digits = compute.bar.diff.digits, format= "f"))]
    comp[, max:= max(dat$max[V1:V2]), .(V1, V2)]
    setorderv(comp, "max")
    comp[, y:= max(dat[.BY, max, on= c("x>=x0", "x<=x1")]), .(x0, x1)]
    # Adjust based on overlaps
    data.table::setorderv(comp, c("y", "x0", "x1"))
    comp[, idx:= .I] # Index to check previous bars
    adj <- strheight("M")*1.2
    for(i in seq(nrow(comp)))
    {
      .c <- sort(comp[comp[i], y, on= c("x0<=x1", "x1>=x0", "idx<=idx")]) # y pos overlapping lines
      .c <- .c[diff(c(.c, Inf))>2*adj  & .c>=comp[i,y]] # y values with enough space
      comp[i, y:= data.table::first(.c)+adj] # Smaller y value + adj
    }
    comp[, {
      segments(x0, y, x1, y, xpd= xpd)
      text(x,
           y,
           var,
           cex= compute.bar.diff.cex,
           offset= 0.1,
           xpd= xpd,
           pos= 3)
    }]
  }

  # Add bar labels ----
  if(!is.null(bar.top.labels))
  {
    text(bar,
         dat$max,
         bar.top.labels,
         pos= 3,
         cex= bar.top.labels.cex,
         xpd= NA)
  }

  # Return bar positions ----
  invisible(bar)
}
