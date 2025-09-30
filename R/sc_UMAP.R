#' Plot UMAP from seurat object
#'
#' Function to plot rasterized umaps from a seurat object.
#'
#' @param umap A seurat object or a data.table/data.frame/matrix containing columns 'cellID', 'umap_1', 'umap_2'.
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

sc_UMAP <- function(umap,
                    value.var= NULL,
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
                    rect.width= NULL,
                    rect.height= NULL,
                    add= FALSE)
{
  # Get UMAP coordinates ----
  coor <- if(class(umap)[1]=="Seurat") {
    as.data.table(umap@reductions$umap@cell.embeddings, keep.rownames = "cellID")
  } else if(is.data.table(umap)) {
    data.table::copy(umap)
  } else
    as.data.table(umap, keep.rownames= "cellID")

  # Checks ----
  if(!all(c("cellID", "umap_1", "umap_2") %in% names(coor)))
    stop("umap should be a Seurat object or a data.table/matrix/data.frame containing columns 'cellID', 'umap_1', 'umap_2'")
  if(class(umap)[1]=="Seurat" && is.null(value.var))
    value.var <- Idents(umap)
  if(class(umap)[1]!="Seurat" && is.null(value.var))
    stop("value.var is missing.")

  # Coerce value.var to numeric OR factor var ----
  # If value.var is a gene symbol
  var <- if(class(umap)[1]=="Seurat" && length(value.var)==1) {
    stopifnot("data" %in% slotNames(umap[[layer]]))
    if(value.var %in% rownames(umap[[layer]]$data))
      umap[[layer]]$data[value.var,] else
        rep(as.numeric(NA), nrow(coor))
  } else if(length(value.var)==nrow(coor) && is.character(value.var)) {
    factor(value.var) # Coerce characters to factors
  } else
    value.var
  if(!(is.factor(var) | is.numeric(var)))
    stop("value.var could not be coerced to numeric or factor.")

  # Make unique data object
  dat <- cbind(coor, data.table(var= var))

  # Select cells ----
  if(!is.null(select.cell))
    dat <- dat[select.cell, on= "cellID"]

  # Default values ----
  # Title
  if(is.null(main))
    main <- if(length(value.var)==1)
      value.var else if(is.numeric(var))
        "Variable" else
          "Group"
  # Color
  if(is.null(col))
    col <- if(is.numeric(var))
      c("skyblue", "red") else
        sc_colors(var, unique.lvl = T)
  # Zlim
  if(is.null(zlim) && is.numeric(var)) {
    z.min <- quantile(var, 0.01, na.rm= T)
    z.max <- quantile(var, 0.99, na.rm= T)
    if(z.max==z.min)
      z.max <- z.min+1
    zlim <- c(z.min, z.max)
    print(zlim)
  }
  # Legend title (not used for factor variable)
  if(is.null(legend.main) && is.numeric(var))
    legend.main <- ifelse(length(value.var)==1, paste(layer,  "counts\n(log2)"), "Variable")

  # Color function ----
  if(is.numeric(var))
    colFUN <- circlize::colorRamp2(zlim, col)

  # Make unique object for plotting ----
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
    ), keyby= .(var, col)]
    if(show.labels) {
      text(labs$cl.x,
           labs$cl.y,
           labs$var,
           cex= cex.lab)
    } else {
      vl_legend(legend = labs$var,
                fill= labs$col,
                cex = cex.legend)
    }
  }

  # Return object
  invisible(coor)
}
