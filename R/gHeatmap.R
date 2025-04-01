#' Create Enhanced Heatmap Visualization
#'
#' @description
#' Creates a highly customizable heatmap using base R graphics, with support
#' for hierarchical clustering, k-means clustering, dendrograms, and cluster
#' visualization. This function provides a comprehensive solution for
#' visualizing matrix data with advanced clustering options.
#'
#' @param x Numeric matrix to be displayed as heatmap.
#' @param cluster.rows Logical or vector specifying row clustering:
#'   * FALSE: no clustering (default)
#'   * TRUE: perform clustering
#'   * vector: pre-defined clustering
#' @param cluster.cols Similar to cluster.rows but for columns.
#' @param kmeans.k Integer specifying number of k-means clusters. NA (default)
#'   uses hierarchical clustering instead.
#' @param zlim Numeric vector of length 2 specifying the minimum and maximum values
#'   for color mapping (e.g., c(-1, 1)). If NULL (default), limits are automatically
#'   calculated based on data range. Note: this argument is ignored if breaks are
#'   specified directly.
#' @param breaks Numeric vector of break points for color mapping. Takes precedence
#'   over `zlim` if both are specified. If NULL (default), breaks are automatically
#'   calculated either from `zlim` if provided, or from the data range. Break points
#'   define the boundaries between colors, so if using n colors, breaks should be a
#'   vector of length n+1.
#' @param col Color palette for the heatmap. If NULL (default), automatically
#'   selects between a diverging blue-white-red palette for data centered around
#'   zero, or a sequential blue-yellow palette for non-negative data.
#' @param show.rownames Logical; whether to show row names (TRUE by default).
#' @param show.colnames Logical; whether to show column names (TRUE by default).
#' @param tilt.colnames Logical; whether to tilt column names (TRUE by default).
#' @param clustering.method Character string specifying hierarchical clustering
#'   method ("complete" by default).
#' @param clustering.distance.rows Character string specifying distance metric
#'   for rows ("euclidean" by default).
#' @param clustering.distance.cols Similar to clustering.distance.rows but for
#'   columns.
#' @param cutree.rows Integer specifying number of row clusters to cut tree into.
#' @param cutree.cols Similar to cutree.rows but for columns.
#' @param show.row.clusters Character specifying position of row cluster
#'   visualization: "right", "left", or FALSE.
#' @param show.col.clusters Character specifying position of column cluster
#'   visualization: "top", "bottom", or FALSE.
#' @param row.clusters.col Vector of two colors for row cluster gradient.
#' @param col.clusters.col Vector of two colors for column cluster gradient.
#' @param legend.cex Numeric scaling factor for legend (1 by default).
#' @param legend.title Character string for legend title ("Value" by default).
#' @param show.numbers Logical or matrix; if TRUE or matrix provided, displays
#'   values in cells.
#' @param numbers.cex Numeric scaling factor for displayed numbers (0.7 by
#'   default).
#' @param cluster.seed Integer seed for reproducible clustering (3453 by
#'   default).
#'
#' @return
#' Invisibly returns a list with two data.tables:
#' * row: contains row names, ordering, and clustering (if applicable)
#' * col: contains column names, ordering, and clustering (if applicable)
#'
#' @details
#' The function provides multiple clustering options:
#' * Hierarchical clustering with various distance metrics
#' * K-means clustering
#' * Pre-defined clustering
#' * Separate clustering for rows and columns
#'
#' Color schemes are automatically selected based on data characteristics:
#' * For data spanning positive and negative: blue-white-red
#' * For non-negative data: blue-yellow
#'
#' Cluster visualization options include:
#' * Dendrograms
#' * Color-coded cluster bars
#' * Cluster labels with sizes
#'
#' @examples
#' # Create example matrix
#' set.seed(1234)
#' test <- matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] <- test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] <- test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] <- test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) <- paste("Test", 1:10, sep = "")
#' rownames(test) <- paste("Gene", 1:20, sep = "")
#'
#'
#' # Basic heatmap
#' gPar(mfrow= c(2,2))
#' gHeatmap(test)
#' title(main= "No clustering")
#'
#' # With hierarchical ordering
#' gHeatmap(test,
#'          cluster.rows= TRUE,
#'          cluster.cols= TRUE,
#'          zlim= c(-5, 5))
#' title(main= "hclust rows & cols")
#' gHeatmap(test,
#'          cluster.rows= TRUE,
#'          cluster.cols= TRUE,
#'          cutree.rows = 3,
#'          cutree.cols = 3,
#'          show.numbers = round(test, 1))
#' title(main= "Cut hclust rows & cols")
#'
#' # Save clusters!
#' cl <- gHeatmap(test2,
#'                cluster.rows= TRUE,
#'                kmeans.k= 5,
#'                cluster.cols= TRUE,
#'                show.row.clusters = "left",
#'                show.col.clusters = "bottom")
#' title(main= "Kmeans rows, hclust cols")
#'
#' # Look at clusters
#' print(cl)
#'
#' @seealso
#' \code{\link{gImage}} for the underlying image plotting
#' \code{\link{gHeatkey}} for the color key
#' \code{\link{hclust}} for hierarchical clustering
#' \code{\link{kmeans}} for k-means clustering
#'
#' @export
gHeatmap <- function(x,
                     cluster.rows= FALSE,
                     cluster.cols= FALSE,
                     kmeans.k= NA,
                     zlim= NULL,
                     breaks= NULL,
                     col= NULL,
                     show.rownames= TRUE,
                     show.colnames= TRUE,
                     tilt.colnames= TRUE,
                     clustering.method = "complete",
                     clustering.distance.rows= "euclidean",
                     clustering.distance.cols= "euclidean",
                     cutree.rows,
                     cutree.cols,
                     show.row.clusters= "right",
                     show.col.clusters= "top",
                     row.clusters.col= c("grey90", "grey40"),
                     col.clusters.col= c("grey90", "grey40"),
                     legend.cex= 1,
                     legend.title= "Value",
                     show.numbers= FALSE,
                     numbers.cex= .7,
                     cluster.seed= 3453)
{
  # Checks and default values ----
  if(is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))
  if(is.null(colnames(x)))
    colnames(x) <- seq(ncol(x))
  if(nrow(x)==1)
    cluster.rows <- FALSE
  if(ncol(x)==1)
    cluster.cols <- FALSE
  if(!show.row.clusters %in% c("left", "right", FALSE))
    stop("show.row.clusters should be one of 'left', 'right' or FALSE.")
  if(!show.col.clusters %in% c("bottom", "top", FALSE))
    stop("show.col.clusters should be one of 'bottom', 'top' or FALSE.")

  # Cluster rows ----
  row.order <- if(isFALSE(cluster.rows)) {
    # No clustering
    seq(nrow(x))
  } else if(length(cluster.rows) == nrow(x)) {
    # User defined clusters
    order(cluster.rows)
  } else if(isTRUE(cluster.rows)) {
    set.seed(cluster.seed)
    if(!is.na(kmeans.k))
    {
      # Kmeans clustering
      cluster.rows <- kmeans(x, centers = kmeans.k)$cluster
      order(cluster.rows)
    } else {
      # Hierachical clustering
      .d <- if(clustering.distance.rows %in% c("pearson", "spearman")) {
        as.dist(1 - cor(t(x),
                        use= "pairwise.complete.obs",
                        method= clustering.distance.rows))
      } else {
        dist(x, method = clustering.distance.rows)
      }
      # Hierarchical clustering
      hcl <- hclust(.d, method = clustering.method)
      # Cutree
      cluster.rows <- if(!missing(cutree.rows)) {
        cutree(hcl, cutree.rows)
      } else
        FALSE
      # Extract dend
      dend <- ggdendro::dendro_data(hcl,
                                    type = "rectangle",
                                    rotate= T)
      rdend <- data.table::as.data.table(dend$segments)
      # Return order
      hcl$order
    }
  } else
    stop("cluster.rows should either match the number of rows in x, or be a logical vector of length 1.")

  # Cluster columns ----
  col.order <- if(isFALSE(cluster.cols)) {
    # No clustering
    seq(ncol(x))
  } else if(length(cluster.cols) == ncol(x)) {
    # User defined clusters
    order(cluster.cols)
  } else if(isTRUE(cluster.cols)) {
    set.seed(cluster.seed)
    .d <- if(clustering.distance.cols %in% c("pearson", "spearman")) {
      as.dist(1 - cor(x,
                      use= "pairwise.complete.obs",
                      method= clustering.distance.cols))
    } else {
      dist(t(x), method = clustering.distance.cols)
    }
    # Hierarchical clustering
    hcl <- hclust(.d, method = clustering.method)
    # Cutree
    cluster.cols <- if(!missing(cutree.cols)) {
      cutree(hcl, cutree.cols)
    } else
      FALSE
    # Extract dend
    dend <- ggdendro::dendro_data(hcl,
                                  type = "rectangle",
                                  rotate= T)
    cdend <- data.table::as.data.table(dend$segments)
    # Return order
    hcl$order
  } else
    stop("cluster.cols should either match the number of cols in x, or be a logical vector of length 1.")

  # Reorder ----
  x <- x[(row.order), , drop= FALSE]
  x <- x[,(col.order), drop= FALSE]
  if(is.matrix(show.numbers)) {
    show.numbers <- show.numbers[(row.order),]
    show.numbers <- show.numbers[,(col.order)]
  }

  # Plot heatmap ----
  obj <- gImage(
    mat = x,
    zlim= zlim,
    breaks = breaks,
    col = col,
    show.rownames = show.rownames && (isFALSE(cluster.rows) || show.row.clusters=="right"),
    show.colnames = show.colnames && (isFALSE(cluster.cols) || show.col.clusters=="top"),
    tilt.colnames = tilt.colnames,
    show.numbers= show.numbers,
    numbers.cex= numbers.cex
  )

  # Border ----
  rect(xleft = par("usr")[1],
       ybottom = par("usr")[3],
       xright = par("usr")[2],
       ytop = par("usr")[4],
       xpd= T)

  # Plot row clusters ----
  # Define line width (will define cluster boxes width)
  line.width <- diff(grconvertX(c(0,1), "line", "user"))
  if(!isFALSE(show.row.clusters) && !isFALSE(cluster.rows)) {
    # Plotting parameters
    rcls <- data.table(cl= cluster.rows[row.order])
    rcls[, name:= paste0(cl, " (n= ", .N, ")"), cl]
    rcls[, top:= .I[1]-0.5, cl]
    rcls[, bottom:= .I[.N]+0.5, cl]
    rcls[, col:= colorRampPalette(row.clusters.col)(.NGRP)[.GRP], keyby= cl]
    setorderv(rcls, "bottom")
    rcls <- unique(rcls)
    # Plot clusters on the right side
    if(show.row.clusters=="right") {
      left <- par("usr")[2]+line.width/5
      rect(xleft = left,
           ybottom = rcls$bottom,
           xright = left+line.width,
           ytop = rcls$top,
           xpd= NA,
           col= rcls$col,
           border= NA)
      text(x = left+line.width/2,
           y = rowMeans(rcls[, .(top, bottom)]),
           labels = rcls$cl,
           cex= par("cex.lab"),
           xpd= NA,
           offset= 0,
           srt= -90)
      # Shift the dendrogram
      left <- left+line.width
    } else if(show.row.clusters=="left") {
      # Or plot on the left side
      axis(2,
           at = rowMeans(rcls[, .(top, bottom)]),
           labels = rcls$name,
           tick = FALSE)
    }
    # Add lines
    abline(h= rcls[, c(top[-1], bottom[-(.N)])])
  }
  # Add row clusters dendrograms ----
  if(exists("rdend")) {
    if(!exists("left"))
      left <- par("usr")[2]
    segments(left+rdend$y/diff(range(rdend[,c(y, yend)]))*line.width,
             par("usr")[4]+rdend$x-0.5,
             left+rdend$yend/diff(range(rdend[,c(y, yend)]))*line.width,
             par("usr")[4]+rdend$xend-0.5,
             xpd= NA)
    left <- left+line.width
  }

  # Plot columns clusters ----
  # Define line height (will define cluster boxes height)
  line.height <- diff(grconvertY(c(0,1), "line", "user"))
  if(!isFALSE(show.col.clusters) && !isFALSE(cluster.cols)) {
    # Plotting parameters
    ccls <- data.table(cl= cluster.cols[col.order])
    ccls[, name:= paste0(cl, "\n(n= ", .N, ")"), cl]
    ccls[, left:= .I[1]-0.5, cl]
    ccls[, right:= .I[.N]+0.5, cl]
    ccls[, col:= colorRampPalette(col.clusters.col)(.NGRP)[.GRP], keyby= cl]
    setorderv(ccls, "left")
    ccls <- unique(ccls)
    # Plot clusters on the top
    if(show.col.clusters=="top") {
      bottom <- par("usr")[4]+line.height/5
      rect(xleft = ccls$left,
           ybottom = bottom,
           xright = ccls$right,
           ytop = bottom+line.height,
           xpd= NA,
           col= ccls$col,
           border= NA)
      text(x = rowMeans(ccls[, .(left, right)]),
           y = bottom+line.height/2,
           cex= par("cex.lab"),
           labels = ccls$cl,
           xpd= NA,
           offset= 0)
      # Shift columns' dendrogram
      bottom <- bottom+line.height
    } else if(show.col.clusters=="bottom") {
      # Or plot on the bottom
      axis(1,
           at = rowMeans(ccls[, .(left, right)]),
           labels = ccls$name,
           tick = FALSE)
    }
    # Add lines
    abline(v= ccls[, c(left[-1], right[-(.N)])])
  }

  # Add columns dendrograms ----
  if(exists("cdend")) {
    if(!exists("bottom"))
      bottom <- par("usr")[4]
    segments(cdend$x,
             bottom+cdend$y/diff(range(cdend[,c(y, yend)]))*line.height,
             cdend$xend,
             bottom+cdend$yend/diff(range(cdend[,c(y, yend)]))*line.height,
             xpd= NA)
  }

  # Add legend
  if(!exists("left"))
    left <- par("usr")[2]
  gHeatkey(col= obj$col,
           breaks = obj$breaks,
           left = left+line.width/2,
           top = par("usr")[4],
           cex = legend.cex,
           height = 4*legend.cex,
           width = 0.75*legend.cex,
           main = legend.title)

  # Return clusters and row orders ----
  row <- data.table(row.name= rownames(x),
                    row.order= row.order)
  if(!isFALSE(cluster.rows))
    row[, row.cluster:= cluster.rows]
  col <- data.table(col.name= colnames(x),
                    col.order= col.order)
  if(!isFALSE(cluster.cols))
    col[, col.cluster:= cluster.cols]

  # Return clusters ----
  obj <- list(row= row, col= col)
  invisible(obj)
}
