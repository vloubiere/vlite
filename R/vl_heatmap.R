#' Create Enhanced Heatmap Visualization
#'
#' Generates a customizable heatmap with support for hierarchical clustering, k-means clustering,
#' dendrograms, and cluster visualization. This function provides advanced options for visualizing
#' matrix data with clustering and color-coded annotations.
#'
#' @param x A matrix, data.table or data.frame of numeric values or factors (that will be coerced to numeric) to be displayed as a heatmap.
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
#' @param row.annotations A vector of length nrow(x) containing row annotations.
#' @param col.annotations  A vector of length ncol(x) containing row annotations.
#' @param show.legend Should the legend be plotted?
#'        Possible values are TRUE, FALSE, "right" (similar to TRUE) or "top". Default= "right".
#' @param legend.title Character string for the legend title. Default= "Value".
#' @param show.numbers Logical or matrix; if set to TRUE or a matrix is provided, displays values in cells.
#' @param show.rownames Logical; whether to display row names. Default= TRUE.
#' @param show.colnames Logical; whether to display column names. Default= TRUE.
#' @param main Add title to the heatmap. Default= NA.
#' @param useRaster logical; if TRUE a bitmap raster is used to plot the image instead of polygons. Default= FALSE.
#' @param clustering.distance.rows Character string specifying the distance metric for rows. Default is "euclidean".
#' @param clustering.distance.cols Similar to clustering.distance.rows but for columns.
#' @param clustering.method Character string specifying the hierarchical clustering method. Default is "complete".
#' @param cutree.rows Integer specifying the number of row clusters to cut the tree into.
#' @param cutree.cols Similar to cutree.rows but for columns.
#' @param show.row.clusters Character specifying the position of row cluster visualization: "right", "left", or FALSE.
#' @param show.col.clusters Character specifying the position of column cluster visualization: "top", "bottom", or FALSE.
#' @param show.row.dendro If rows are clustered using hclust, should the dendrogram be shown? Default= TRUE.
#' @param show.col.dendro If cols are clustered using hclust, should the dendrogram be shown? Default= TRUE.
#' @param gap.width The width of the gap between clusters, expressed as a fraction of the plot limits. Default= 1/40.
#' @param cluster.seed Integer seed for reproducible clustering. Default= 3453.
#' @param row.clusters.col A vector of two colors for the row cluster gradient. Default is c("grey90", "grey40").
#' @param col.clusters.col A vector of two colors for the column cluster gradient. Default is c("grey90", "grey40").
#' @param row.annotations.col Row annotations colors. Default= rainbow(9)[1:7].
#' @param col.annotations.col Col annotations colors. Default= rainbow(9)[1:7].
#' @param na.col Color for NA values. Default= "ligthgrey".
#' @param legend.cex Numeric scaling factor for the legend. Default= 1.
#' @param numbers.cex Numeric scaling factor for the size of displayed numbers. Default= 0.7.
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
                       row.annotations= NULL,
                       col.annotations= NULL,
                       show.legend= "right",
                       legend.title= "Value",
                       show.numbers= FALSE,
                       show.rownames= TRUE,
                       show.colnames= TRUE,
                       main= NA,
                       useRaster= FALSE,
                       clustering.distance.rows= "euclidean",
                       clustering.distance.cols= "euclidean",
                       clustering.method = "complete",
                       cutree.rows,
                       cutree.cols,
                       show.row.clusters= "right",
                       show.col.clusters= "top",
                       show.row.dendro= TRUE,
                       show.col.dendro= TRUE,
                       gap.width= 1/40,
                       cluster.seed= 3453,
                       row.clusters.col= c("grey90", "grey40"),
                       col.clusters.col= c("grey90", "grey40"),
                       row.annotations.col= rainbow(9)[1:7],
                       col.annotations.col= rainbow(9)[1:7],
                       na.col= "lightgrey",
                       legend.cex= 1,
                       numbers.cex= .7)
{
  # Checks and default values ----
  # Coerce x to numeric matrix (useful for factors)
  if(!is.matrix(x)) {
    checkClass <- unique(sapply(x, class))
    if(length(checkClass)>1)
      stop("x cannot contain mixed classes.")
    if(checkClass=="factor") {
      allLvls <- unique(c(sapply(x, levels)))
      # Factors to numeric
      x <- lapply(x, function(x) as.numeric(factor(x, allLvls)))
      x <- do.call(cbind, lapply(x, as.numeric))
    }
    # Numeric matrix
    x <- as.matrix(x)
  } else {
    checkClass <- "numeric"
  }
  if(is.logical(x))
    x <- apply(x, 2, as.numeric)
  if(!is.numeric(x))
    stop("x should contain numeric values, factors or logical values.")
  if(is.null(colnames(x)))
    colnames(x) <- seq(ncol(x))
  if(is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))
  if(nrow(x)==1)
    cluster.rows <- FALSE
  if(ncol(x)==1)
    cluster.cols <- FALSE
  if(!is.null(row.annotations) && length(row.annotations) != nrow(x))
    stop("row.annotations should be a vector of length nrow(x)")
  if(!is.null(row.annotations) && !is.factor(row.annotations))
    row.annotations <- factor(row.annotations, unique(row.annotations))
  if(!is.null(col.annotations) && length(col.annotations) != ncol(x))
    stop("col.annotations should be a vector of length ncol(x)")
  if(!is.null(col.annotations) && !is.factor(col.annotations))
    col.annotations <- factor(col.annotations, unique(col.annotations))
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
  if(!isFALSE(cluster.rows) && show.row.clusters=="left")
    show.rownames <- FALSE
  if(!isFALSE(cluster.cols) && show.col.clusters=="bottom")
    show.colnames <- FALSE
  if(isTRUE(show.legend))
    show.legend <- "right"
  if(isTRUE(show.numbers))
    show.numbers <- x
  if(!isFALSE(show.numbers) & !is.matrix(show.numbers))
    show.numbers <- as.matrix(show.numbers)

  # Default breaks ----
  if(is.null(breaks)) {
    breaks <- if(checkClass=="factor") {
      # Factor
      seq(min(x, na.rm= TRUE),
          max(x, na.rm= TRUE))
    } else if(min(x, na.rm= TRUE)<0 & max(x, na.rm= TRUE)>0) {
      lims <- max(abs(x), na.rm= TRUE)
      # Centered on 0
      seq(-lims,
          lims,
          length.out= 21)
    } else {
      # Unique sign
      seq(min(x, na.rm= TRUE),
          max(x, na.rm= TRUE),
          length.out= 21)
    }
  }

  # Default colors ----
  if(is.null(col)) {
    col <- if(checkClass=="factor") {
      # Factor
      rainbow(length(allLvls))
    } else if(breaks[1] < 0 & breaks[length(breaks)] > 0) {
      # Positive and negative values
      c("royalblue1", "white", "red")
    } else {
      # Unique sign
      c("blue", "yellow")
    }
  }
  col <- colorRampPalette(col)(length(breaks))

  # Compute row and column gap widths ----
  if(nrow(x)>=ncol(x)) {
    row.gap.width <- nrow(x)*gap.width
    col.gap.width <- row.gap.width*(ncol(x)/nrow(x))
  }else{
    col.gap.width <- ncol(x)*gap.width
    row.gap.width <- col.gap.width*(nrow(x)/ncol(x))
  }

  # Cluster rows ----
  rows <- heatmap.get.clusters(name= rownames(x),
                               idx= seq(nrow(x)),
                               annot = row.annotations,
                               x = x,
                               clusters = cluster.rows,
                               distance = clustering.distance.rows,
                               method = clustering.method,
                               cutree = cutree.rows,
                               kmeans.k = kmeans.k,
                               cluster.seed = cluster.seed,
                               cluster.col = row.clusters.col,
                               annot.col = row.annotations.col,
                               gap.width = row.gap.width)
  rdend <- rows$dend
  rows <- rows$obj
  setnames(rows,
           c("idx", "im", "pos"),
           c("line.idx", "im.col", "y.pos"))

  # Cluster columns ----
  cols <- heatmap.get.clusters(name= colnames(x),
                               idx= seq(ncol(x)),
                               annot = col.annotations,
                               x = t(x),
                               clusters = cluster.cols,
                               distance = clustering.distance.cols,
                               method = clustering.method,
                               cutree = cutree.cols,
                               kmeans.k = NA,
                               cluster.seed = cluster.seed,
                               cluster.col = col.clusters.col,
                               annot.col = col.annotations.col,
                               gap.width = col.gap.width)
  cdend <- cols$dend
  cols <- cols$obj
  setnames(cols,
           c("idx", "im", "pos"),
           c("column.idx", "im.row", "x.pos"))

  # Reorder matrix based on clustering ----
  x <- x[(rows$order), , drop=FALSE]
  x <- x[, (cols$order), drop=FALSE]
  if(is.matrix(show.numbers)) {
    show.numbers <- show.numbers[(rows$order), , drop=FALSE]
    show.numbers <- show.numbers[, (cols$order), drop= FALSE]
  }

  # Prepare image for plotting ----
  # Transpose
  im <- t(x)
  # Clip outlier values to min/max color breaks
  im[im<min(breaks)] <- min(breaks, na.rm= TRUE)
  im[im>max(breaks)] <- max(breaks, na.rm= TRUE)
  # Set NA values to min(breaks)-1
  im[is.na(im)] <- min(breaks)-1
  # Adjust breaks to afford NAs
  NAbreaks <- c(breaks[1]-c(2,1), breaks)
  NAcols <- c(na.col, col)

  # Order rows and cols before plotting ----
  rows <- rows[(order)]
  cols <- cols[(order)]

  # Initiate full image ----
  image(x= 1,
        y= 1,
        z= matrix(NA),
        breaks= c(-1,1),
        col= NA,
        xlim= cols[c(1, .N), x.pos]+c(-0.5, 0.5),
        ylim= rows[c(.N, 1), y.pos]+c(0.5, -0.5),
        xlab= NA,
        ylab= NA,
        axes= FALSE)

  # Plot clusters iteratively ----
  rows[, {
    cols[, {
      image(x= if(length(x.pos)==1) x.pos+c(-0.5, 0.5) else x.pos,
            y= if(length(y.pos)==1) y.pos+c(-0.5, 0.5) else y.pos,
            z= im[im.row, im.col, drop= FALSE],
            breaks= NAbreaks,
            col= NAcols,
            xlab= NA,
            ylab= NA,
            axes= FALSE,
            useRaster= useRaster,
            add= TRUE)
      rect(xleft = x.pos[1]-0.5,
           ybottom = rev(y.pos)[1]+0.5,
           xright = rev(x.pos)[1]+0.5,
           ytop = y.pos[1]-0.5,
           xpd= NA)
    }, cluster]
  }, cluster]

  # Plot numbers ----
  if(!isFALSE(show.numbers)) {
    text(x= rep(cols$x.pos, each= nrow(rows)),
         y= rep(rows$y.pos, nrow(cols)),
         labels = c(show.numbers),
         cex= numbers.cex)
  }

  # Plot axes ----
  # X axis
  if(show.rownames) {
    axis(2,
         rows$y.pos,
         rows$name,
         lwd = 0,
         lwd.ticks = 0)
  }
  # Y axis
  if(show.colnames) {
    # Check if labels will overlap
    cn.width <- strwidth(colnames(x))/2
    cn.n <- seq(ncol(x))
    colnames.ov <- (cn.n-cn.width)[-1] < (cn.n+cn.width)[-length(cn.n)]
    # If yes, tilt them
    if(any(colnames.ov)) {
      tiltAxis(x= cols$x.pos,
               labels= cols$name)
    } else {
      # Otherwise, default axis
      axis(1,
           at = cols$x.pos,
           labels = cols$name,
           lwd = 0,
           lwd.ticks = 0,
           padj= -1.25)
    }
  }

  # Compute margins positions and define line width/height (used as ref) ----
  right.mar <- par("usr")[2]
  top.mar <- par("usr")[4]
  line.width <- diff(grconvertX(c(0,1), "line", "user"))
  line.height <- diff(grconvertY(c(0,1), "line", "user"))

  # Plot row annotations ----
  if(!is.null(row.annotations)) {
    # Adjust margin
    right.mar <- right.mar+line.width/5
    # Plot
    rows[, {
      rect(xleft = right.mar,
           ybottom = rev(y.pos)[1]+0.5,
           xright = right.mar+line.width,
           ytop = y.pos[1]-0.5,
           col= annot.col[1],
           xpd= NA,
           border= NA)
    }, rleid(annot, cluster)]
    # Adjust margin
    right.mar <- right.mar+line.width
  }

  # Plot col annotations ----
  if(!is.null(col.annotations)) {
    # Adjust margin
    top.mar <- top.mar+line.height/5
    # Plot
    cols[, {
      rect(xleft = x.pos[1]-0.5,
           ybottom = top.mar,
           xright = rev(x.pos)[1]+0.5,
           ytop = top.mar+line.height,
           col= annot.col[1],
           xpd= NA,
           border= NA)
    }, rleid(annot, cluster)]
    # Adjust margin
    top.mar <- top.mar+line.height
  }

  # Plot row clusters ----
  if(sum(!is.na(rows$cluster)) && !isFALSE(show.row.clusters)) {
    # Plot on the right side
    if(show.row.clusters=="right") {
      # Adjust margin
      right.mar <- right.mar+line.width/5
      # Plot
      rows[, {
        rect(xleft = right.mar,
             ybottom = rev(y.pos)[1]+0.5,
             xright = right.mar+line.width,
             ytop = y.pos[1]-0.5,
             col= cluster.col[1],
             xpd= NA,
             border= NA)
        text(x = right.mar+line.width/2,
             y = mean(y.pos),
             labels = cluster[1],
             cex= par("cex.lab"),
             xpd= NA,
             offset= 0,
             srt= -90)
      }, rleid(cluster)]
      # Adjust margin
      right.mar <- right.mar+line.width
    } else if(show.row.clusters=="left") {
      # Plot on the left side
      axis(2,
           at = rows[, mean(y.pos), cluster]$V1,
           labels = rows[, paste0(cluster, " (n= ", formatC(.N, big.mark = ","), ")"), cluster]$V1,
           tick = FALSE)
    }
  }

  # Plot column clusters ----
  if(sum(!is.na(cols$cluster)) && !isFALSE(show.col.clusters)) {
    # On the top
    if(show.col.clusters=="top") {
      # Adjust margin
      top.mar <- top.mar+line.height/5
      # Plot
      cols[, {
        rect(xleft = x.pos[1]-0.5,
             ybottom = top.mar,
             xright = rev(x.pos)[1]+0.5,
             ytop = top.mar+line.height,
             col= cluster.col[1],
             xpd= NA,
             border= NA)
        text(x = mean(x.pos),
             y = top.mar+line.height/2,
             labels = cluster[1],
             cex= par("cex.lab"),
             xpd= NA,
             offset= 0)
      }, rleid(cluster)]
      # Adjust margin
      top.mar <- top.mar+line.height
    } else if(show.col.clusters=="bottom") {
      # Or plot on the bottom
      axis(1,
           at = cols[, mean(x.pos), cluster]$V1,
           labels = cols[, paste0(cluster, "\nn= ", formatC(.N, big.mark = ",")), cluster]$V1,
           tick = FALSE)
    }
  }

  # Add dendrograms ----
  # Rows
  if(exists("rdend") && show.row.dendro) {
    # For each segment, interpolate new y0 and y1
    rdend[, y0 := {
      start <- rows$y.pos[floor(x)]
      end   <- rows$y.pos[ceiling(x)]
      start + (end - start) * (x - floor(x))
    }]
    rdend[, y1 := {
      start <- rows$y.pos[floor(xend)]
      end   <- rows$y.pos[ceiling(xend)]
      start + (end - start) * (xend - floor(xend))
    }]
    # Now plot
    segments(
      right.mar + rdend$y    / diff(range(rdend[,c(y, yend)])) * line.width,
      par("usr")[4] + rdend$y0    - 0.5,
      right.mar + rdend$yend / diff(range(rdend[,c(y, yend)])) * line.width,
      par("usr")[4] + rdend$y1    - 0.5,
      xpd = NA
    )
    # Adjust margin
    right.mar <- right.mar + line.width
  }

  # Columns
  if(exists("cdend") && show.col.dendro) {
    # For each segment, interpolate new x0 and x1
    cdend[, x0 := {
      start <- cols$x.pos[floor(x)]
      end   <- cols$x.pos[ceiling(x)]
      start + (end - start) * (x - floor(x))
    }]
    cdend[, x1 := {
      start <- cols$x.pos[floor(xend)]
      end   <- cols$x.pos[ceiling(xend)]
      start + (end - start) * (xend - floor(xend))
    }]
    # Now plot
    segments(cdend$x0,
             top.mar+cdend$y/diff(range(cdend[,c(y, yend)]))*line.height,
             cdend$x1,
             top.mar+cdend$yend/diff(range(cdend[,c(y, yend)]))*line.height,
             xpd= NA)
    # Adjust margin
    top.mar <- top.mar+line.height
  }

  # Add heatkey ----
  if(!isFALSE(show.legend)) {
    # Adjust plotting position
    adj.x <- ifelse(show.legend=="top",
                    0,
                    (right.mar-par("usr")[2])/line.width-.5)
    adj.y <- ifelse(show.legend=="top",
                    (top.mar-par("usr")[4])/line.height-.5,
                    0)
    # Main heatmap heatkey
    heatkey(col= col,
            breaks = if(checkClass=="factor") factor(allLvls, allLvls) else breaks,
            position = show.legend,
            adj.x = adj.x,
            adj.y = adj.y,
            cex = legend.cex,
            main = legend.title)
    # Adjust top margin (for title)
    top.mar <- top.mar+((show.legend=="top")*3*line.height)
    # Row annotations
    if(!is.null(row.annotations)) {
      # Adjust pos
      adj.y <- adj.y-ifelse(show.legend=="right", 5.5, -2.5)*legend.cex
      # Add heatkey
      rann <- unique(rows[,.(annot, annot.col)])
      setorderv(rann, "annot")
      heatkey(col= rann$annot.col,
              breaks = factor(rann$annot, unique(rann$annot)),
              position = show.legend,
              adj.x = adj.x,
              adj.y = adj.y,
              cex = legend.cex,
              main = "Rows")
    }
    # Col annotations
    if(!is.null(col.annotations)) {
      # Adjust plotting position
      adj.y <- adj.y-ifelse(show.legend=="right", 5.5, -2.5)*legend.cex
      # Add heatkey
      cann <- unique(cols[,.(annot, annot.col)])
      setorderv(cann, "annot")
      heatkey(col= cann$annot.col,
              breaks = factor(cann$annot, unique(cann$annot)),
              position = show.legend,
              adj.x = adj.x,
              adj.y = adj.y,
              cex = legend.cex,
              main = "Columns")
    }
  }

  # Add title ----
  if(!is.na(main))
    title(main= main,
          line = max(c(1, (top.mar-par("usr")[4])/line.height+.25)))

  # Return clusters ----
  obj <- list(rows= rows[order(line.idx), .(name, line.idx, cluster, order, y.pos)],
              cols= cols[order(column.idx), .(name, column.idx, cluster, order, x.pos)])
  invisible(obj)
}
