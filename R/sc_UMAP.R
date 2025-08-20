#' Plot UMAP from seurat object
#'
#' Function to plot rasterized umaps from a seurat object.
#'
#' @param dat A seurat object.
#' @param value.var The value to plot. Should either be a factor vector matching the number of cells in the object,
#' or a unique gene symbol.
#' @param show.labels Should the factor levels be shown on the UMAP itself? Otherwise, they will be shown in a legend.
#' @param main
#' @param col
#' @param NA.col
#' @param zlim
#' @param legend.main
#' @param cex.pts
#' @param transparency
#' @param cex.lab Expansion factor for labels.
#' @param cex.legend Expansion factor for the legend.
#' @param randomize
#' @param layer The layer to use. Default= "SCT".
#'
#' @return Invisibly returns a list containing the ordered heatmap matrix.
#' @export

sc_UMAP <- function(dat,
                    value.var= Idents(dat),
                    show.labels= TRUE,
                    main= NULL,
                    col= NULL,
                    NA.col= "lightgrey",
                    zlim= NULL,
                    legend.main= NULL,
                    cex.pts= .4,
                    pch= 16,
                    transparency= .7,
                    cex.lab= .7,
                    cex.legend= .7,
                    randomize= FALSE,
                    select.cell= NULL,
                    layer= "SCT",
                    add= FALSE)
{
  # Checks ----
  if(class(dat)!="Seurat")
    stop("dat should be a seurat object.")
  if(!length(value.var) %in% c(1, ncol(dat[[layer]])))
    stop("value.var should either match the number of cells in the object or be a unique gene symbol.")
  if(length(value.var)==1 & !"counts" %in% slotNames(dat[[layer]]))
    stop(paste0("The 'counts' slot is missing from the ", layer, " layer."))

  # Determine type of plot ----
  plot.type <- fcase(length(value.var)==1, "Gene",
                     is.numeric(value.var), "Variable",
                     default = "Group")

  # Retrieve layer ----
  SCT <- dat[[layer]]$counts

  # Extract cell UMAP coordinates ----
  coor <- as.data.table(dat@reductions$umap@cell.embeddings, keep.rownames = T)

  # Select cells ----
  if(!is.null(select.cell)) {
    SCT <- SCT[,select.cell]
    value.var <- droplevels(value.var[select.cell])
    coor <- coor[select.cell, ]
  }

  # Retrieve variable ----
  var <- if(length(value.var)==ncol(SCT)){
    if(is.character(value.var))
      factor(value.var, unique(value.var)) else
        value.var
  } else {
    value.var <- as.character(value.var)
    if(value.var %in% rownames(SCT))
      SCT[value.var,] else
        rep(as.numeric(NA), ncol(SCT))
  }

  # Defaults ----
  if(is.null(main))
    main <- ifelse(plot.type=="Gene", value.var, plot.type)
  if(is.null(col))
    col <- if(is.numeric(var))
      c("skyblue", "red") else
        sc_colors(var, unique.lvl = T)
  if(is.null(zlim) && is.numeric(var)) {
    z.min <- min(c(0, quantile(var, 0.01, na.rm= T)), na.rm= T)
    z.max <- max(c(1, quantile(var, 0.99, na.rm= T)), na.rm= T)
    if(z.max==z.min)
      z.max <- z.min+1
    zlim <- c(z.min, z.max)
  }
  if(is.null(legend.main) && is.numeric(var))
    legend.main <- ifelse(plot.type=="Gene", paste(layer,  "counts\n(log2)"), "Variable")

  # Compute color function ----
  if(is.numeric(var))
    colFUN <- circlize::colorRamp2(zlim, col)

  # Add variables to coordinates ----
  coor[, var:= var]
  if(is.numeric(var))
    coor[, col:= ifelse(is.na(var), NA.col, colFUN(as.numeric(var)))] else
      coor[, col:= col[as.numeric(var)]]

  # Order ----
  if(randomize) {
    coor <- coor[sample(.N)]
  } else if(is.numeric(var))
    setorderv(coor, "var")

  # Plot ----
  if(add) {
    points(
      coor$umap_1,
      coor$umap_2,
      cex= cex.pts,
      col = adjustcolor(coor$col, transparency),
      pch= pch
    )
  } else {
    rasterScatterplot(
      coor$umap_1,
      coor$umap_2,
      cex= cex.pts,
      col = adjustcolor(coor$col, transparency),
      xlab= "UMAP 1",
      ylab= "UMAP 2",
      pch= pch
    )
  }

  # Title
  title(main= main)

  # Labels
  if(is.numeric(var)) {
    br <- seq(zlim[1], zlim[2], length.out = 21)
    heatkey(breaks = br,
            col = colFUN(br),
            main = legend.main,
            cex = cex.legend)
  } else {
    labs <- coor[, .(
      cl.x= median(umap_1),
      cl.y= median(umap_2)
    ), .(var, col)]
    if(show.labels) {
      text(labs$cl.x,
           labs$cl.y,
           labs$var,
           cex= cex.lab)
    } else {
      vl_legend("topright",
                legend = labs$var,
                fill= labs$col,
                cex = cex.legend)
    }
  }

  # Return object
  invisible(coor)
}
