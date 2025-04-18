#' Create Enhanced Heatmap Visualization
#'
#' Generates a customizable heatmap with support for hierarchical clustering, k-means clustering,
#' dendrograms, and cluster visualization. This function provides advanced options for visualizing
#' matrix data with clustering and color-coded annotations.
#'
#' @param x A numeric matrix to be displayed as a heatmap.
#' @param cluster.rows Logical or vector specifying row clustering:
#'   * TRUE: perform clustering (default)
#'   * FALSE: no clustering
#'   * vector: pre-defined clustering
#' @param cluster.cols Similar to cluster.rows but for columns. Default= FALSE.
#' @param kmeans.k Integer specifying the number of k-means clusters. Defaults= NA (uses hierarchical clustering).
#' @param breaks A numeric vector specifying the breakpoints for color mapping.
#'   Defaults to 21 evenly spaced breaks spanning the range of x.
#' @param col A vector of colors corresponding to the values in breaks. If set to NULLs,
#'   a diverging blue-white-red palette is used for data spanning positive and negative values,
#'   or a sequential blue-yellow palette.
#' @param show.rownames Logical; whether to display row names. Default= TRUE.
#' @param show.colnames Logical; whether to display column names. Default= TRUE.
#' @param tilt.colnames Logical; whether to tilt column names for better readability. Default= TRUE.
#' @param main Ad title for the heatmap. Default= NA.
#' @param clustering.method Character string specifying the hierarchical clustering method. Default is "complete".
#' @param clustering.distance.rows Character string specifying the distance metric for rows. Default is "euclidean".
#' @param clustering.distance.cols Similar to clustering.distance.rows but for columns.
#' @param cutree.rows Integer specifying the number of row clusters to cut the tree into.
#' @param cutree.cols Similar to cutree.rows but for columns.
#' @param show.row.clusters Character specifying the position of row cluster visualization: "right", "left", or FALSE.
#' @param show.col.clusters Character specifying the position of column cluster visualization: "top", "bottom", or FALSE.
#' @param show.row.dendro If rows are clustered using hclust, should the dendrogram be shown? Default= TRUE.
#' @param show.col.dendro If cols are clustered using hclust, should the dendrogram be shown? Default= TRUE.
#' @param show.legend Should the legend be plotted? Possible values are TRUE, FALSE, "right" (similar to TRUE) or "top". Default= "right".
#' @param row.clusters.col A vector of two colors for the row cluster gradient. Default is c("grey90", "grey40").
#' @param col.clusters.col A vector of two colors for the column cluster gradient. Default is c("grey90", "grey40").
#' @param na.col Color for NA values. Default= "ligthgrey".
#' @param legend.cex Numeric scaling factor for the legend. Default= 1.
#' @param legend.title Character string for the legend title. Default= "Value".
#' @param show.numbers Logical or matrix; if set to TRUE or a matrix is provided, displays values in cells.
#' @param numbers.cex Numeric scaling factor for the size of displayed numbers. Default= 0.7.
#' @param cluster.seed Integer seed for reproducible clustering. Default= 3453.
#' @param useRaster logical; if TRUE a bitmap raster is used to plot the image instead of polygons. Default= FALSE.
#'
#' @return
#' Invisibly returns a list with two `data.table` objects:
#' * row: Contains row names, ordering, and clustering (if applicable).
#' * col: Contains column names, ordering, and clustering (if applicable).
#'
#' @details
#' The function supports:
#' - Hierarchical clustering with customizable distance metrics and methods.
#' - K-means clustering for rows or columns.
#' - Pre-defined clustering for rows or columns.
#' - Visualization of dendrograms and color-coded cluster bars.
#'
#' Color schemes are automatically selected based on the data:
#' - For data spanning positive and negative values: blue-white-red.
#' - For non-negative data: blue-yellow.
#'
#' The default behavior for `breaks` is to create 21 evenly spaced intervals spanning the range of `x`.
#' If `x` contains both positive and negative values, the range is symmetrically extended around zero.
#'
#' @examples
#' # Create example matrix
#' set.seed(1234)
#' mat <- matrix(rnorm(200), 20, 10)
#' colnames(mat) <- paste("Sample", 1:10, sep = "")
#' rownames(mat) <- paste("Gene", 1:20, sep = "")
#'
#' # Basic heatmap
#' vl_par()
#' vl_heatmap(mat)
#'
#' # Hierarchical clustering
#' vl_heatmap(mat, cluster.rows = TRUE, cluster.cols = TRUE)
#'
#' # K-means clustering
#' vl_heatmap(mat, kmeans.k = 3, cluster.cols = TRUE)
#'
#' # Display cell values
#' vl_heatmap(mat, show.numbers = round(mat, 1))
#'
#' @seealso
#' \code{\link{vl_image}} for the underlying image plotting.
#' \code{\link{hclust}} for hierarchical clustering.
#' \code{\link{kmeans}} for k-means clustering.
#'
#' @export
vl_heatmap <- function(x,
                       cluster.rows= TRUE,
                       cluster.cols= FALSE,
                       kmeans.k= NA,
                       breaks= NULL,
                       col= NULL,
                       show.rownames= TRUE,
                       show.colnames= TRUE,
                       tilt.colnames= TRUE,
                       main= NA,
                       clustering.method = "complete",
                       clustering.distance.rows= "euclidean",
                       clustering.distance.cols= "euclidean",
                       cutree.rows,
                       cutree.cols,
                       show.row.clusters= "right",
                       show.col.clusters= "top",
                       row.clusters.col= c("grey90", "grey40"),
                       col.clusters.col= c("grey90", "grey40"),
                       na.col= "lightgrey",
                       show.row.dendro= TRUE,
                       show.col.dendro= TRUE,
                       show.legend= "right",
                       legend.cex= 1,
                       legend.title= "Value",
                       show.numbers= FALSE,
                       numbers.cex= .7,
                       cluster.seed= 3453,
                       useRaster= FALSE)
{
  # Checks and default values ----
  if(!is.matrix(x))
    x <- as.matrix(x)
  if(is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))
  if(is.null(colnames(x)))
    colnames(x) <- seq(ncol(x))
  if(nrow(x)==1)
    cluster.rows <- FALSE
  if(ncol(x)==1)
    cluster.cols <- FALSE
  if(!show.row.clusters %in% c(TRUE, FALSE, "left", "right"))
    stop("show.row.clusters should be one of TRUE, FALSE, 'left' or 'right'.")
  if(isTRUE(show.row.clusters))
    show.row.clusters <- "left"
  if(!show.col.clusters %in% c(TRUE, FALSE, "bottom", "top"))
    stop("show.col.clusters should be one of TRUE, FALSE, 'bottom' or 'top'.")
  if(isTRUE(show.col.clusters))
    show.col.clusters <- "top"
  if(!show.legend %in% c(TRUE, FALSE, "right", "top"))
    stop("show.col.clusters should be one of TRUE, FALSE, 'right' or 'top'.")
  if(isTRUE(show.legend))
    show.legend <- "right"

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
  obj <- vl_image(
    mat = x,
    breaks = breaks,
    col = col,
    show.rownames = show.rownames && (isFALSE(cluster.rows) || show.row.clusters=="right"),
    show.colnames = show.colnames && (isFALSE(cluster.cols) || show.col.clusters=="top"),
    tilt.colnames = tilt.colnames,
    show.numbers= show.numbers,
    numbers.cex= numbers.cex,
    useRaster= useRaster,
    na.col = na.col
  )

  # Border ----
  rect(xleft = par("usr")[1],
       ybottom = par("usr")[3],
       xright = par("usr")[2],
       ytop = par("usr")[4],
       xpd= T)

  # Compute margsins positions ----
  right.mar <- par("usr")[2]
  top.mar <- par("usr")[4]

  # Define line width and height, which will serve as ref below ----
  line.width <- diff(grconvertX(c(0,1), "line", "user"))
  line.height <- diff(grconvertY(c(0,1), "line", "user"))

  # Plot row clusters ----
  if(!isFALSE(cluster.rows)) {

    # Plotting parameters
    rcls <- data.table(cl= cluster.rows[row.order])
    rcls[, name:= paste0(cl, " (n= ", .N, ")"), cl]
    rcls[, top:= .I[1]-0.5, cl]
    rcls[, bottom:= .I[.N]+0.5, cl]
    rcls[, col:= colorRampPalette(row.clusters.col)(.NGRP)[.GRP], keyby= cl]
    setorderv(rcls, "bottom")
    rcls <- unique(rcls)

    # Add cluster lines
    abline(h= rcls[, c(top[-1], bottom[-(.N)])])

    # Add cluster labels
    if(!isFALSE(show.row.clusters)) {

      # Plot clusters on the right side
      if(show.row.clusters=="right") {
        right.mar <- right.mar+line.width/5
        rect(xleft = right.mar,
             ybottom = rcls$bottom,
             xright = right.mar+line.width,
             ytop = rcls$top,
             xpd= NA,
             col= rcls$col,
             border= NA)
        text(x = right.mar+line.width/2,
             y = rowMeans(rcls[, .(top, bottom)]),
             labels = rcls$cl,
             cex= par("cex.lab"),
             xpd= NA,
             offset= 0,
             srt= -90)

        # Shift margin
        right.mar <- right.mar+line.width

      } else if(show.row.clusters=="left") {

        # Or plot on the left side
        axis(2,
             at = rowMeans(rcls[, .(top, bottom)]),
             labels = rcls$name,
             tick = FALSE)
      }
    }
  }

  # Plot columns clusters ----
  if(!isFALSE(cluster.cols)) {

    # Plotting parameters
    ccls <- data.table(cl= cluster.cols[col.order])
    ccls[, name:= paste0(cl, "\n(n= ", .N, ")"), cl]
    ccls[, left:= .I[1]-0.5, cl]
    ccls[, right:= .I[.N]+0.5, cl]
    ccls[, col:= colorRampPalette(col.clusters.col)(.NGRP)[.GRP], keyby= cl]
    setorderv(ccls, "left")
    ccls <- unique(ccls)

    # Add cluster lines
    abline(v= ccls[, c(left[-1], right[-(.N)])])

    if(!isFALSE(show.col.clusters)) {

      # Plot clusters on the top
      if(show.col.clusters=="top") {
        top.mar <- top.mar+line.height/5
        rect(xleft = ccls$left,
             ybottom = top.mar,
             xright = ccls$right,
             ytop = top.mar+line.height,
             xpd= NA,
             col= ccls$col,
             border= NA)
        text(x = rowMeans(ccls[, .(left, right)]),
             y = top.mar+line.height/2,
             cex= par("cex.lab"),
             labels = ccls$cl,
             xpd= NA,
             offset= 0)

        # Shift margin
        top.mar <- top.mar+line.height

      } else if(show.col.clusters=="bottom") {

        # Or plot on the bottom
        axis(1,
             at = rowMeans(ccls[, .(left, right)]),
             labels = ccls$name,
             tick = FALSE)
      }
    }
  }

  # Add dendrograms ----
  # Rows
  if(exists("rdend") && show.row.dendro) {
    segments(right.mar+rdend$y/diff(range(rdend[,c(y, yend)]))*line.width,
             par("usr")[4]+rdend$x-0.5,
             right.mar+rdend$yend/diff(range(rdend[,c(y, yend)]))*line.width,
             par("usr")[4]+rdend$xend-0.5,
             xpd= NA)

    # Shift margin
    right.mar <- right.mar+line.width
  }
  # Columns
  if(exists("cdend") && show.col.dendro) {
    segments(cdend$x,
             top.mar+cdend$y/diff(range(cdend[,c(y, yend)]))*line.height,
             cdend$xend,
             top.mar+cdend$yend/diff(range(cdend[,c(y, yend)]))*line.height,
             xpd= NA)

    # Shift margin
    top.mar <- top.mar+line.height
  }

  # Add legend ----
  if(!isFALSE(show.legend)) {
    heatkey(col= obj$col,
            breaks = obj$breaks,
            position = show.legend,
            adj.x = ifelse(show.legend=="top",
                           0,
                           (right.mar-par("usr")[2])/line.width-.5),
            adj.y = ifelse(show.legend=="top",
                           (top.mar-par("usr")[4])/line.height-.5,
                           0),
            cex = legend.cex,
            main = legend.title)
    top.mar <- top.mar+((show.legend=="top")*3*line.height)
  }


  # Add title ----
  title(main= main,
        line = (top.mar-par("usr")[4])/line.height+.5)

  # Return clusters and row orders ----
  row <- data.table(row.name= rownames(x)[order(row.order)],
                    row.order= row.order)
  if(!isFALSE(cluster.rows))
    row[, row.cluster:= cluster.rows]
  col <- data.table(col.name= colnames(x)[order(col.order)],
                    col.order= col.order)
  if(!isFALSE(cluster.cols))
    col[, col.cluster:= cluster.cols]

  # Return clusters ----
  obj <- list(row= row, col= col)
  invisible(obj)
}
