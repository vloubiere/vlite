#' Add labels heatmap
#'
#' Add axis with repelled labels to a plot.
#'
#' @param side The axis is placed as follows: 2 = left (default), 4 = right.
#' @param at Y position of the labels.
#' @param labels Labels
#' @param cex cex expansion factor for the axis. Default= par("cex.axis").
#' @param width width of the labels. Default= 1.
#'
#' @returns
#' @export
#'
#' @examples
addRepelAxis <- function(
    side= 2,
    at,
    labels,
    cex= par("cex.axis"),
    width= 1
) 
{
  stopifnot(length(at)==length(labels))
  stopifnot(side %in% c(2,4))
  width <- width*2
  
  # Make data table ----
  dat <- data.table(at= at, lab= labels)
  dat <- dat[!is.na(at) & !is.na(lab)]
  dat[, row_id := .I]
  dat <- dat[order(at)]
  
  # Compute grid of positions ----
  all.dist <- diff(par("usr")[c(3,4)])
  all.dist <- floor(all.dist/(strheight("(M", cex= cex)*1.5))
  all.dist <- seq(par("usr")[3], par("usr")[4], length.out= all.dist)
  if(length(all.dist)<nrow(dat))
    all.dist <- seq(par("usr")[3], par("usr")[4], length.out= nrow(dat))
  all.dist <- sort(all.dist) # Sort to prevent line crossings
  
  # GPT script >>>>>>>>>>>>>>>
  # gives the globally optimal non-crossing assignment, minimizing total distance while preserving order.
  # Dimensions
  n <- nrow(dat)
  m <- length(all.dist)
  # Distance between each label anchor and each possible grid position
  cost <- outer(dat$at, all.dist, function(a, g) abs(a - g))
  # Dynamic programming tables
  dp <- matrix(Inf, n, m)
  prev <- matrix(NA_integer_, n, m)
  # First label can use any grid position
  dp[1, ] <- cost[1, ]
  # Find optimal ordered assignment
  for(i in 2:n) {
    best_prev <- Inf
    best_j <- NA_integer_
    
    for(j in 2:m) {
      if(dp[i - 1, j - 1] < best_prev) {
        best_prev <- dp[i - 1, j - 1]
        best_j <- j - 1
      }
      
      dp[i, j] <- cost[i, j] + best_prev
      prev[i, j] <- best_j
    }
  }
  # Recover chosen grid positions
  idx <- integer(n)
  idx[n] <- which.min(dp[n, ])
  for(i in n:2) {
    idx[i - 1] <- prev[i, idx[i]]
  }
  # Store adjusted positions
  dat[, adj.at := all.dist[idx]]
  # Restore original order
  dat <- dat[order(row_id)]
  # GPT script <<<<<<<<<<<<<<<<
  
  # When the original position has no overlap, keep it
  dat[, text.start:= adj.at+strheight("(M", cex= cex)*0.75]
  dat[, text.end:= adj.at-strheight("(M", cex= cex)*0.75]
  dat[, box.start:= min(c(text.start, text.end, at)), .(text.start, text.end)]
  dat[, box.end:= max(c(text.start, text.end, at)), .(text.start, text.end)]
  dat[, orig.start:= at+strheight("(M", cex= cex)*0.75]
  dat[, orig.end:= at-strheight("(M", cex= cex)*0.75]
  dat[, o.start:= min(c(orig.start, orig.end)), .(orig.start, orig.end)]
  dat[, o.end:= max(c(orig.start, orig.end)), .(orig.start, orig.end)]
  check <- dat[dat, .N, .EACHI, on= c("box.end>=o.start", "box.start<=o.end")]$N
  dat[check==1, adj.at:= at]
  
  # Axis plotting coordinates ----
  pos <- if(side==4)
    par("usr")[2] else if(side==2)
      par("usr")[1]
  marg <- par("mgp")[2]*diff(grconvertX(c(0,1), from = "line", "user"))*width
  
  # Plot labels ----
  text(
    x = pos+ifelse(side==4, marg, -marg),
    y = dat$adj.at,
    labels = dat$lab,
    pos= side,
    xpd= NA,
    offset= 0.25,
    cex= cex
  )
  
  # Plot lines ----
  dat[, {
    lines(
      c(
        pos,
        pos+ifelse(side==4, marg, -marg)/4,
        pos+ifelse(side==4, marg, -marg)*3/4,
        pos+ifelse(side==4, marg, -marg)
      ),
      c(
        at,
        at,
        adj.at,
        adj.at
      ),
      xpd= NA
    )
  }, .(at, adj.at)]
}