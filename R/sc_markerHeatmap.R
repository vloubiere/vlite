#' Plot and Cluster Top Marker Genes as a Heatmap from Single-Cell Transcriptome Data
#'
#' Selects top marker genes per cluster based on multiple statistical cutoffs, then visualizes their expression as a heatmap.
#' Optionally clusters rows and/or orders columns, and allows highlighting of selected genes.
#'
#' @param dat A data.table containing marker gene statistics (see ?sc_computeMarkerGenes).
#' @param topN Number of top genes to select per cluster. Ddefault= 3.
#' If set to 0 or NULL, only the genes listed in selection (next argument) are used.
#' @param selection Character vector of selected gene symbols that should be added to the heatmap.
#' Default= NULL.
#' @param cluster.selection Should select genes be used to cluster rows or not? Default= FALSE.
#' @param order.selection Should selected gene columns be ordered?
#' @param cluster.rows Logical; whether to cluster heatmap rows. Default= TRUE.
#' @param cutree.rows Number of row clusters to cut the dendrogram into. Default= 10.
#' @param show.col.clusters Should column cluster names be shown? Default= TRUE.
#' @param perc.cutoff Minimum percentage of expressing cells within the cluster. Default= 10.
#' @param perc.diff.cutoff Minimum difference in percentage of expressing cells within vs. outside the cluster.
#' Default= 20.
#' @param order.var Variable to order genes by; one of "perc.diff", "auc", "logFC", or "log2OR".
#' Default= "perc.diff".
#' @param padj.wilcox.cutoff Adjusted Wilcoxon p-value cutoff. Default= 0.05.
#' @param logFC.cutoff Log fold change cutoff. Default= 0.
#' @param padj.fisher.cutoff Adjusted Fisher p-value cutoff. Default= NULL.
#' @param log2OR.cutoff Log2 odds ratio cutoff. Default= NULL.
#' @param main Title.
#' @param auc.cutoff AUC cutoff. Default= NULL.
#' @param alternative "greater" for upregulated markers, "smaller" for downregulated. Default= "greater".
#' @param add.cells.number If set to TRUE, a barplot specifying the number of cells is plotted. Default= TRUE.
#' @param pdf.file If specified, the heatmap will be saved in a pdf with an optimized layout. Default= NULL.
#'
#' @return Invisibly returns a list containing the ordered heatmap matrix.
#' @export
sc_markerHeatmap <- function(
    dat,
    topN= 3,
    selection= NULL,
    cluster.selection= FALSE,
    order.selection= FALSE,
    cluster.rows= FALSE,
    cutree.rows= NULL,
    show.col.clusters= TRUE,
    row.annotation.column= NULL,
    row.annotation.colors= c("#A3D8FF","#B6FF8A", "#FFE066", "#FF6F69"),
    row.annotation.title= NULL,
    col.annotation.column= NULL,
    col.annotation.colors= c("lightgrey", "chartreuse", "red", "gold", "royalblue1", "royalblue3", "blue"),
    col.annotation.title= NULL,
    perc.cutoff= 10,
    perc.diff.cutoff= 20,
    order.var= "perc.diff",
    padj.wilcox.cutoff= 0.05,
    logFC.cutoff= 0.25,
    padj.fisher.cutoff= NULL,
    log2OR.cutoff= NULL,
    auc.cutoff= NULL,
    main= NA,
    alternative= "greater",
    add.cells.number= TRUE,
    pdf.file= NULL
) {
  # Checks ----
  if(is.null(topN))
    topN <- 0
  if(!"Nclust" %in% names(dat))
    stop("dat should contain a Nclust column (number of cells within cluster).")
  if(!is.data.table(dat))
    stop("dat should be a data.table containing marker genes.")
  if(!is.null(col.annotation.column) && !col.annotation.column %in% names(dat))
    stop("col.annotation.column missing from data.table.")
  if(!is.null(row.annotation.column) && !row.annotation.column %in% names(dat))
    stop("row.annotation.column missing from data.table.")
  if(!is.null(row.annotation.column) && is.null(row.annotation.title))
    row.annotation.title <- row.annotation.column
  if(!is.null(col.annotation.column) && is.null(col.annotation.title))
    col.annotation.title <- col.annotation.column
  # Check if selected genes exist
  if(!is.null(selection)) {
    selection <- unique(selection)
    # Check missing genes
    if(!all(selection %in% dat$symbol)) {
      warning(
        paste("The following genes are missing from the table:",
              paste0(setdiff(selection, dat$symbol), collapse = ", "))
      )
    }
    selection <- intersect(selection, dat$symbol)
    if(length(selection)==0)
      selection <- NULL
  }
  if(is.null(selection))
    selection <- character()

  # Extract top markers ----
  top <- sc_topMarkers(
    dat= dat,
    topN= topN,
    perc.cutoff= perc.cutoff,
    perc.diff.cutoff= perc.diff.cutoff,
    padj.wilcox.cutoff= padj.wilcox.cutoff,
    logFC.cutoff= logFC.cutoff,
    padj.fisher.cutoff= padj.fisher.cutoff,
    log2OR.cutoff= log2OR.cutoff,
    auc.cutoff= auc.cutoff,
    order.var= order.var,
    alternative= alternative
  )

  # Retrieve top cluster of selected genes ----
  sel <- dat[selection, on= "symbol"]
  sel <- sel[, .SD[which.max(perc)], symbol]

  # Save selection names (add a star before plotting labels) ----
  selection.names <- sel$symbol

  # Check ----
  if(nrow(sel)+nrow(top)==0)
    stop("No genes selected. Either relax filtering paramters or add a selection of genes.")

  # If selected genes should be used for clustering ----
  if(cluster.selection) {
    # Remove genes that were already found as top markers
    sel <- sel[!symbol %in% top$symbol]
    # Add them to the list of top markers
    top <- rbind(top, sel, fill= TRUE)
    # Empty selection
    sel <- sel[0]
  }

  # Retrieve number of cells ----
  Ncells <- dat[, Nclust[1], keyby= cluster]$V1

  # Retrieve row annotation ----
  if(!is.null(row.annotation.column))
    row.annotation.column <- dat[, get(row.annotation.column)[1], keyby= cluster]$V1

  # If top/selected genes are selected ----
  if(nrow(top)) {
    # Get % of expressing cells for top marker genes ----
    mat <- dcast(dat[symbol %in% top$symbol],
                 cluster~symbol,
                 value.var = "perc",
                 drop= TRUE)

    # Casting
    mat <- as.matrix(mat, 1)

    # Cluster rows ----
    if(isTRUE(cluster.rows)) {

      # Clustering
      cl <- vl_heatmap(mat,
                       cluster.rows = cluster.rows,
                       clustering.distance.rows = "spearman",
                       plot = F)

      # Reorder objects
      mat <- mat[cl$rows$order,]
      top[, cluster:= droplevels(factor(cluster, unique(rownames(mat))))]
      Ncells <- Ncells[cl$rows$order]
      if(!is.null(row.annotation.column))
        row.annotation.column <- row.annotation.column[cl$rows$order]
    }

    # Order columns ----
    setorderv(top, "cluster", 1)
    mat <- mat[, unique(top$symbol), drop= F]
  }

  # If selected genes should be added to the matrix ----
  if(nrow(sel)) {

    # Get percentage of expressing cells
    mat.sel <- dcast(dat[symbol %in% selection],
                     cluster~symbol,
                     value.var = "perc",
                     drop= TRUE)

    # Casting
    mat.sel <- as.matrix(mat.sel, 1)

    # Order rows based on clustering
    if(nrow(top) && isTRUE(cluster.rows))
      mat.sel <- mat.sel[cl$rows$order,]

    # Order columns
    mat.sel <- if(order.selection) {
      sel[, cluster:= droplevels(factor(cluster, rownames(mat.sel)))]
      setorderv(sel, "cluster", 1)
      mat.sel[, unique(sel$symbol), drop= F]
    } else
      mat.sel[, sel$symbol, drop= F]
  }

  # Get column clusters ----
  ccl <- droplevels(top[, cluster[1], symbol, drop= F]$V1)
  if(nrow(sel)) {
    selection.cl.name <- rev(make.unique(c(levels(ccl), "selection")))[1]
    levels(ccl) <- c(levels(ccl), selection.cl.name)
    ccl <- c(ccl, factor(rep(selection.cl.name, ncol(mat.sel))))
  }

  # Make final matrix for heatmap ----
  hm.mat <- if(nrow(top) & nrow(sel))
    cbind(mat, mat.sel) else if(nrow(top))
      mat else if(nrow(sel))
        mat.sel

  # Retrieve columns annotations ----
  if(!is.null(col.annotation.column)) {
    if(is.null(col.annotation.title))
      col.annotation.title <- col.annotation.column
    idx <- match(colnames(hm.mat), dat$symbol)
    col.annotation.column <- dat[[col.annotation.column]][idx]
    if(!is.factor(col.annotation.column))
      col.annotation.column <- factor(col.annotation.column, unique(col.annotation.column))
  }

  # Add an * to selected gene names ----
  old.colnames <- colnames(hm.mat)
  idx <- colnames(hm.mat) %in% selection.names
  colnames(hm.mat)[idx] <- paste0(colnames(hm.mat)[idx], "*")

  # Heatmap ----
  cl <- vl_heatmap(
    hm.mat,
    cluster.rows = if(nrow(top)==0) FALSE else if(isTRUE(cluster.rows)) mat else cluster.rows,
    clustering.distance.rows = "spearman",
    cutree.rows = cutree.rows,
    cluster.cols = ccl,
    show.col.clusters = show.col.clusters,
    row.annotations = row.annotation.column,
    row.annotations.col = adjustcolor(row.annotation.colors, .7),
    row.annotations.title = row.annotation.title,
    col.annotations = col.annotation.column,
    col.annotations.col = adjustcolor(col.annotation.colors, .7),
    col.annotations.title = col.annotation.title,
    col= c("white", "red"),
    legend.title = "% exp. cells",
    col.gap.width = .5,
    row.gap.width = .5,
    legend.cex = .6,
    show.numbers = round(hm.mat),
    numbers.cex= .4,
    pdf.file = pdf.file,
    pdf.close = !add.cells.number,
    main= main
  )

  # Add barplot ----
  if(add.cells.number) {
    lab.width <- strwidth(paste0("    ", rownames(hm.mat)), cex= par("cex.axis"))
    right <- par("usr")[1]-max(lab.width)
    line.width <- diff(grconvertX(c(0,1), "line", "user"))
    norm.width <- Ncells/max(Ncells)*line.width
    rect(xleft = right-norm.width,
         ybottom = cl$rows$y.pos-.4,
         xright= right,
         ytop = cl$rows$y.pos+.4,
         xpd= NA,
         border= NA,
         col= "lightgrey")
    text(right-norm.width,
         cl$rows$y.pos,
         formatC(Ncells, big.mark = ","),
         pos= 2,
         cex= par("cex.axis"),
         offset= 0.2,
         xpd= NA)
    if(!is.null(pdf.file))
      dev.off()
  }

  # Return matrix ----
  colnames(hm.mat) <- old.colnames
  hm.mat <- hm.mat[cl$rows$order, cl$cols$order]
  invisible(hm.mat)
}
