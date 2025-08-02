#' Create Enhanced Heatmap Visualization
#'
#' Generates a customizable heatmap with support for hierarchical clustering, k-means clustering,
#' dendrograms, and cluster visualization. This function provides advanced options for visualizing
#' matrix data with clustering and color-coded annotations.
#'
#' @param x A matrix, data.table or data.frame of numeric values or factors (that will be coerced to numeric) to be displayed as a heatmap.
#' @param cluster.rows Logical or vector specifying row clustering:
#'   * TRUE: perform clustering (default).
#'   * FALSE: no clustering.
#'   * vector: pre-defined clustering. Will be coerced to factors.
#'   * data.table/matrix: clustering will be performed on this subset instead of x.
#' @param cluster.cols Similar to cluster.rows but for columns (of note, 'matrix' option is not available here).
#' Default= FALSE.
#' @param kmeans.k Integer specifying the number of k-means clusters for rows. Defaults= NA (uses hierarchical clustering).
#' @param breaks A numeric vector specifying the breakpoints for color mapping.
#'   Can either be the length of col (centered breaks) or one element longer (in which case breaks will correspond to edges).
#'   If set to NULL, defaults to the different levels for factors (centered breaks), a diverging scale
#'   of 21 values (centered on 0) or to 20 values spanning the range of x (breaks as edges).
#' @param col A vector of colors corresponding to the values in breaks.
#'   Can either be the length of breaks (centered breaks) or one element shorter (breaks used as edges).
#'   If set to NULL, rainbow colors are used for factors, a diverging blue-white-red palette for data
#'   spanning positive and negative values, or a sequential blue-yellow palette for data with a unique sign.
#' @param row.annotations A vector of length nrow(x) containing row annotations.
#' @param col.annotations  A vector of length ncol(x) containing row annotations.
#' @param show.legend Should the legend be plotted?
#'        Possible values are TRUE, FALSE, "right" (similar to TRUE) or "top". Default= "right".
#' @param legend.title Character string for the legend title. Default= "Value".
#' @param row.annotations.title Title for the row annotations legend. Default= "Rows".
#' @param col.annotations.title Title for the column annotations legend. Default= "Columns".
#' @param show.numbers Logical or matrix; if set to TRUE or a matrix is provided, displays values in cells.
#' @param show.rownames Logical; whether to display row names. Default= TRUE.
#' @param show.colnames Logical; whether to display column names. Default= TRUE.
#' @param main Add title to the heatmap. Default= NA.
#' @param useRaster logical; if TRUE a bitmap raster is used to plot the image instead of polygons. Default= FALSE.
#' @param clustering.distance.rows Character string specifying the distance metric for rows. Default is "euclidean".
#' @param clustering.distance.cols Similar to clustering.distance.rows but for columns.
#' @param clustering.method Character string specifying the hierarchical clustering method. Default is "complete".
#' @param cutree.rows Integer specifying the number of row clusters to cut the tree into. Default= NULL.
#' @param cutree.cols Similar to cutree.rows but for columns. Default= NULL.
#' @param show.row.clusters Character specifying the position of row cluster visualization: "right", "left", or FALSE.
#' @param show.col.clusters Character specifying the position of column cluster visualization: "top", "bottom", or FALSE.
#' @param show.row.dendro If rows are clustered using hclust, should the dendrogram be shown? Default= TRUE.
#' @param show.col.dendro If cols are clustered using hclust, should the dendrogram be shown? Default= TRUE.
#' @param row.gap.width The height of the gap between row clusters, expressed as a fraction of rows' height. Default= 1/10.
#' @param col.gap.width The width of the gap between column clusters, expressed as a fraction of columns' height. Default= 1/10.
#' @param cluster.seed Integer seed for reproducible clustering. Default= 3453.
#' @param row.clusters.col A vector of two colors for the row cluster gradient. Default is c("grey90", "grey40").
#' @param col.clusters.col A vector of two colors for the column cluster gradient. Default is c("grey90", "grey40").
#' @param row.annotations.col Row annotations colors. Default= rainbow(9)[1:7].
#' @param col.annotations.col Col annotations colors. Default= rainbow(9)[1:7].
#' @param na.col Color for NA values. Default= "ligthgrey".
#' @param legend.cex Numeric scaling factor for the legend. Default= 1.
#' @param numbers.cex Numeric scaling factor for the size of displayed numbers. Default= 0.7.
#' @param show.grid Should the matric grid be shown? Default= FALSE.
#' @param grid.lwd Line width used for the grid.
#' @param pdf.file If specified, the heatmap will be saved into the provided pdf file, with an optimized and nice layout.
#' Default= NULL
#' @param pdf.par If TRUE, the plotting parametersof the PDF will be automatically set to nice default values. Oterwise,
#' the current par() settings will be used. Default= TRUE.
#'
#' @return
#' Invisibly returns a list with two `data.table` objects:
#' * row: Contains row names, ordering, and clustering (if applicable).
#' * col: Contains column names, ordering, and clustering (if applicable).
#'
#' @examples
#' # Create example matrix
#' set.seed(1234)
#' mat <- matrix(rnorm(200), 20, 10)
#' colnames(mat) <- paste("Sample", 1:10, sep = "")
#' rownames(mat) <- paste("Gene", 1:20, sep = "")
#' mat <- round(mat, 2)
#'
#' # Basic heatmap
#' vl_par()
#' vl_heatmap(mat)
#'
#' # Hierarchical clustering
#' vl_heatmap(mat,
#' cluster.rows = T,
#' cluster.cols = T,
#' cutree.rows = 3,
#' cutree.cols = 3,
#' show.numbers = T)
#'
#' @export
vl_heatmap <- function(x,
                       cluster.rows= TRUE,
                       cluster.cols= FALSE,
                       kmeans.k= NA,
                       kmeans.order= FALSE,
                       breaks= NULL,
                       col= NULL,
                       row.annotations= NULL,
                       col.annotations= NULL,
                       show.legend= "right",
                       legend.title= "Value",
                       row.annotations.title= "Rows",
                       col.annotations.title= "Columns",
                       show.numbers= FALSE,
                       show.rownames= TRUE,
                       show.colnames= TRUE,
                       main= NA,
                       useRaster= FALSE,
                       clustering.distance.rows= "euclidean",
                       clustering.distance.cols= "euclidean",
                       clustering.method = "complete",
                       cutree.rows= NULL,
                       cutree.cols= NULL,
                       show.row.clusters= "right",
                       show.col.clusters= "top",
                       show.row.dendro= TRUE,
                       show.col.dendro= TRUE,
                       row.gap.width= 1/10,
                       col.gap.width= 1/10,
                       cluster.seed= 3453,
                       row.clusters.col= c("grey90", "grey40"),
                       col.clusters.col= c("grey90", "grey40"),
                       row.annotations.col= rainbow(9)[1:7],
                       col.annotations.col= rainbow(9)[1:7],
                       na.col= "lightgrey",
                       legend.cex= 1,
                       numbers.cex= .7,
                       show.grid= FALSE,
                       grid.lwd= .25,
                       plot= T,
                       pdf.file= NULL,
                       pdf.par= TRUE)
{
  # Coerce x to numeric matrix (useful for factors) ----
  x <- toNumMatrix(x)
  checkClass <- x$checkClass
  x <- x$x
  if(is.null(colnames(x)))
    colnames(x) <- seq(ncol(x))
  if(is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))

  # Check clustering parameters values ----
  # Rows
  if(nrow(x)==1 && isTRUE(cluster.rows))
    cluster.rows <- FALSE
  if(!isTRUE(cluster.rows) && !isFALSE(cluster.rows)) {
    if(is.vector(cluster.rows))
      cluster.rows <- factor(cluster.rows)
    if(is.factor(cluster.rows)) {
      if(length(cluster.rows) != nrow(x))
        stop("cluster.rows should match the number of rows in x.")
    } else if(is.matrix(cluster.rows)) {
      cx <- toNumMatrix(cluster.rows)$x
      if(nrow(cx)!=nrow(x))
        stop("cluster.rows matrix should have the same number of rows as x.")
      rownames(cx) <- rownames(x)
      cluster.rows <- TRUE
    } else
      stop("cluster.rows format not recognized.")
  }
  if(!is.logical(kmeans.order) || kmeans.order>1)
    stop("kmeans.order should be TRUE or FALSE.")
  if(isTRUE(show.row.clusters))
    show.row.clusters <- "left"
  if(!show.row.clusters %in% c(FALSE, "left", "right"))
    stop("show.row.clusters should be one of TRUE, FALSE, 'left' or 'right'.")
  if(isTRUE(cluster.rows) && clustering.distance.rows %in% c("pearson", "spearman") && ncol(x)==1)
    stop("Error: Rows cannot be clustered using pearson or spearman when ncol(x)==1")
  # Columns
  if(ncol(x)==1 && isTRUE(cluster.cols))
    cluster.cols <- FALSE
  if(!isTRUE(cluster.cols) && !isFALSE(cluster.cols)) {
    if(is.vector(cluster.cols))
      cluster.cols <- factor(cluster.cols)
    if(is.factor(cluster.cols)) {
      if(length(cluster.cols) != ncol(x))
        stop("cluster.cols should match the number of columns in x.")
    } else
      stop("cluster.cols format not recognized")
  }
  if(isTRUE(show.col.clusters))
    show.col.clusters <- "top"
  if(!show.col.clusters %in% c(FALSE, "bottom", "top"))
    stop("show.col.clusters should be one of TRUE, FALSE, 'bottom' or 'top'.")
  if(isTRUE(cluster.cols) && clustering.distance.cols %in% c("pearson", "spearman") && nrow(x)==1)
    stop("Error: Columns cannot be clustered using pearson or spearman when nrow(x)==1")

  # Check annotations ----
  if(!is.null(row.annotations)){
    if(is.vector(row.annotations))
      row.annotations <- factor(row.annotations)
    if(is.factor(row.annotations)) {
      if(length(row.annotations) != nrow(x))
        stop("row.annotations should match the number of rows in x.")
    } else
      stop("row.annotations format not recognized.")
  }
  if(!is.null(col.annotations)){
    if(is.vector(col.annotations))
      col.annotations <- factor(col.annotations)
    if(is.factor(col.annotations)) {
      if(length(col.annotations) != ncol(x))
        stop("col.annotations should match the number of columns in x.")
    } else
      stop("col.annotations format not recognized.")
  }

  # Check legend parameters ----
  if(isTRUE(show.legend))
    show.legend <- "right"
  if(!show.legend %in% c(FALSE, "right", "top"))
    stop("show.col.clusters should be one of TRUE, FALSE, 'right' or 'top'.")
  if(!isFALSE(cluster.rows) && show.row.clusters=="left")
    show.rownames <- FALSE
  if(!isFALSE(cluster.cols) && show.col.clusters=="bottom")
    show.colnames <- FALSE
  if(isTRUE(show.numbers))
    show.numbers <- x
  if(!isFALSE(show.numbers) & !is.matrix(show.numbers))
    show.numbers <- as.matrix(show.numbers)

  # Default breaks ----
  if(is.null(breaks)) {
    breaks <- if(checkClass=="factor") {
      # Factors (centered breaks)
      seq(max(x, na.rm= TRUE))
    } else if(min(x, na.rm= TRUE)<0 & max(x, na.rm= TRUE)>0) {
      # Centered on 0 (centered breaks)
      lims <- max(abs(x), na.rm= TRUE)
      seq(-lims,
          lims,
          length.out= 21)
    } else {
      # Unique sign (breaks as edges)
      seq(min(x, na.rm= TRUE),
          max(x, na.rm= TRUE),
          length.out= 20)
    }
  }

  # Default colors ----
  if(is.null(col)) {
    col <- if(checkClass=="factor") {
      # Factor (centered breaks)
      colorRampPalette(rainbow(9)[1:7])(length(breaks))
    } else if(breaks[1] < 0 & breaks[length(breaks)] > 0) {
      # Positive and negative values (centered breaks)
      colorRampPalette(c("royalblue1", "white", "red"))(length(breaks))
    } else {
      # Unique sign (breaks as edges)
      colorRampPalette(c("blue", "yellow"))(length(breaks)-1)
    }
  }

  # Center breaks if they are the same length as col ----
  if(length(breaks)==length(col)) {
    d <- diff(breaks)
    breaks <- c(
      breaks[1] - d[1]/2,
      (breaks[-1] + breaks[-length(breaks)])/2,
      breaks[length(breaks)] + d[length(d)]/2
    )
  } else {
    # By default, breaks are used as edges
    col <- colorRampPalette(col)(length(breaks)-1)
  }

  # Cluster columns object ----
  cols <- heatmap.get.clusters(dim= "col",
                               x = x,
                               annot = col.annotations,
                               clusters = cluster.cols,
                               distance = clustering.distance.cols,
                               method = clustering.method,
                               cutree = cutree.cols,
                               kmeans.k = NA, # Kmeans not implemented for columns
                               cluster.seed = cluster.seed,
                               cluster.col = col.clusters.col,
                               annot.col = col.annotations.col,
                               gap.width = col.gap.width)
  cdend <- cols$dend
  cols <- cols$obj

  # Cluster rows object ----
  rows <- heatmap.get.clusters(dim= "row",
                               x = if(exists("cx")) cx else x,
                               annot = row.annotations,
                               clusters = cluster.rows,
                               distance = clustering.distance.rows,
                               method = clustering.method,
                               cutree = cutree.rows,
                               kmeans.k = kmeans.k,
                               kmeans.cl.names = if(kmeans.order) cols$name else NULL,
                               cluster.seed = cluster.seed,
                               cluster.col = row.clusters.col,
                               annot.col = row.annotations.col,
                               gap.width = row.gap.width)
  rdend <- rows$dend
  rows <- rows$obj

  # Clip outlier values to min/max color breaks ----
  x[x<min(breaks)] <- min(breaks, na.rm= TRUE)
  x[x>max(breaks)] <- max(breaks, na.rm= TRUE)
  # Set NA values to min(breaks)-1
  x[is.na(x)] <- min(breaks)-1
  # Adjust breaks to afford NAs
  NAbreaks <- c(breaks[1]-c(2,1), breaks[-1])
  NAcols <- c(na.col, col)

  # Initiate pdf file ----
  if(!is.null(pdf.file)) {
    # Compute adjusted number of rows heatmap
    Nrows <- nrow(x)
    if(!isFALSE(cluster.rows)) {
      Nrows <- Nrows+(uniqueN(rows$cluster)-1)*row.gap.width
    }
    # Compute adjusted number of columns heatmap
    Ncols <- ncol(x)
    if(!isFALSE(cluster.cols)) {
      Ncols <- Ncols+(uniqueN(cols$cluster)-1)*col.gap.width
    }
    # Compute adjusted number of rows heatmap
    pdf(pdf.file,
        width= Ncols*0.12+4,
        height = Nrows*0.12+4)
    if(pdf.par) {
      old.par <- par()
      vl_par(mai= c(2,2,2,2))
    }
  }

  # Plotting ----
  if(plot) {
    # Plot heatmap ----
    image(x= 1,
          y= 1,
          z= matrix(NA),
          breaks= c(-1,1),
          col= NA,
          xlim= range(cols$x.pos)+c(-0.5, 0.5),
          ylim= rev(range(rows$y.pos))+c(0.5, -0.5), # Firs line at y= 1...
          xlab= NA,
          ylab= NA,
          axes= FALSE)

    # Plot clusters iteratively ----
    rows[, {
      cols[, {
        image(x= if(length(x.pos)==1) x.pos+c(-0.5, 0.5) else x.pos,
              y= if(length(y.pos)==1) y.pos+c(-0.5, 0.5) else y.pos,
              z= t(x[line.idx, column.idx, drop= FALSE]),
              breaks= NAbreaks,
              col= NAcols,
              xlab= NA,
              ylab= NA,
              axes= FALSE,
              useRaster= useRaster,
              add= TRUE)
        if(show.grid) {
          x <- seq(x.pos[1]-0.5, rev(x.pos)[1]+0.5)
          segments(x,
                   rev(y.pos)[1]+0.5,
                   x,
                   y.pos[1]-0.5,
                   lwd= grid.lwd)
          y <- seq(y.pos[1]-0.5, rev(y.pos)[1]+0.5)
          segments(x.pos[1]-0.5,
                   y,
                   rev(x.pos)[1]+0.5,
                   y,
                   lwd= grid.lwd)
        }
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
           labels = c(show.numbers[rows$line.idx, cols$column.idx, drop= FALSE]),
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
      cn.width <- strwidth(cols$name)/2
      colnames.ov <- (cols$x.pos-cn.width)[-1] < (cols$x.pos+cn.width)[-length(cn.width)]
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

    # Compute margins positions and define line width/height (used for plotting) ----
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
             ybottom = y.pos[.N]+0.5,
             xright = right.mar+line.width*.5,
             ytop = y.pos[1]-0.5,
             col= annot.col[1],
             xpd= NA,
             border= NA)
      }, rleid(annot, cluster)]
      # Adjust margin
      right.mar <- right.mar+line.width*.5
    }

    # Plot col annotations ----
    if(!is.null(col.annotations)) {
      # Adjust margin
      top.mar <- top.mar+line.height/5
      # Plot
      cols[, {
        rect(xleft = x.pos[1]-0.5,
             ybottom = top.mar,
             xright = x.pos[.N]+0.5,
             ytop = top.mar+line.height*.5,
             col= annot.col[1],
             xpd= NA,
             border= NA)
      }, rleid(annot, cluster)]
      # Adjust margin
      top.mar <- top.mar+line.height*.5
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
               ybottom = y.pos[.N]+0.5,
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
        }, cluster]
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
               xright = x.pos[.N]+0.5,
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
        }, cluster]
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

    # Add rows dendrogram ----
    if(!is.null(rdend) && show.row.dendro) {
      segments(right.mar + rdend$x0 * line.width,
               par("usr")[4] + rdend$y0,
               right.mar + rdend$x1 * line.width,
               par("usr")[4] + rdend$y1,
               xpd = NA,
               lend= 2)
      # Adjust margin
      right.mar <- right.mar + line.width
    }

    # Add columns dendrogram ----
    if(!is.null(cdend) && show.col.dendro) {
      segments(cdend$x0,
               top.mar+cdend$y0*line.height,
               cdend$x1,
               top.mar+cdend$y1*line.height,
               xpd= NA,
               lend= 2)
      # Adjust margin
      top.mar <- top.mar+line.height
    }

    # Add heatkeys ----
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
              breaks = breaks,
              labels = if(checkClass=="factor") allLvls else NULL,
              position = show.legend,
              adj.x = adj.x,
              adj.y = adj.y,
              cex = legend.cex,
              main = legend.title)
      # Adjust top margin
      top.mar <- top.mar+((show.legend=="top")*3*line.height)

      # Row annotations
      if(!is.null(row.annotations)) {
        # Adjust plotting position
        adj.y <- adj.y-ifelse(show.legend=="right", 5.5, -2.5)*legend.cex
        # Annotations heatkey
        rann <- factor(levels(rows$annot), levels(rows$annot))
        heatkey(breaks = -as.numeric(rann),
                col= row.annotations.col[rann],
                labels = levels(rann),
                position = show.legend,
                adj.x = adj.x,
                adj.y = adj.y,
                cex = legend.cex,
                main = row.annotations.title)
      }

      # Col annotations
      if(!is.null(col.annotations)) {
        # Adjust plotting position
        adj.y <- adj.y-ifelse(show.legend=="right", 5.5, -2.5)*legend.cex
        # Annotations heatkey
        cann <- factor(levels(cols$annot), levels(cols$annot))
        heatkey(breaks = -as.numeric(cann),
                col= col.annotations.col[cann],
                labels = levels(cann),
                position = show.legend,
                adj.x = adj.x,
                adj.y = adj.y,
                cex = legend.cex,
                main = col.annotations.title)
      }
    }

    # Add title ----
    if(!is.na(main))
      title(main= main,
            line = max(c(1, (top.mar-par("usr")[4])/line.height+.25)))
  }

  # Close pdf ----
  if(!is.null(pdf.file)) {
    dev.off()
    print(paste("PDF file ->", pdf.file))
  }

  # Return clusters ----
  obj <- list(rows= rows[order(line.idx), .(name, line.idx, cluster, order, y.pos)],
              cols= cols[order(column.idx), .(name, column.idx, cluster, order, x.pos)])
  invisible(obj)
}
