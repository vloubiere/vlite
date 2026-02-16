hc_polygon_fixed <- function(object, ir = NULL, x1 = NULL, x2 = NULL,
                             gp = grid::gpar(),
                             end_type = c("expanding", "average", "shrinking")) {
  end_type <- match.arg(end_type)[1]
  
  if (object@MODE == "pixel") {
    oi <- HilbertCurve:::.ENV$I_PLOT
    grid::seekViewport(paste0("hilbert_curve_", HilbertCurve:::.ENV$I_PLOT))
    
    hc2 <- HilbertCurve::HilbertCurve(
      s = object@data_range[1],
      e = object@data_range[2],
      mode = "normal",
      level = min(object@LEVEL, 9),
      newpage = FALSE,
      zoom = object@ZOOM,
      legend = FALSE,  # FIX: only once
      start_from = object@start_from,
      first_seg = object@first_seg,
      padding = grid::unit(0, "mm")
    )
    
    # call the fixed function again so it works for both pixel/non-pixel paths
    hc_polygon_fixed(hc2, ir = ir, x1 = x1, x2 = x2, gp = gp, end_type = end_type)
    
    grid::seekViewport(name = paste0("hilbert_curve_", oi, "_global"))
    grid::upViewport()
    return(invisible(NULL))
  }
  
  grid::seekViewport(name = paste0("hilbert_curve_", HilbertCurve:::get_plot_index()))
  
  polygons <- HilbertCurve:::get_polygons(
    object, ir = ir, x1 = x1, x2 = x2, end_type = end_type
  )
  
  gp <- HilbertCurve:::validate_gpar(
    gp, default = list(lty = 1, lwd = 1, col = 1, fill = "transparent")
  )
  
  if (length(gp$lty) == 1)  gp$lty  <- rep(gp$lty,  length(polygons))
  if (length(gp$lwd) == 1)  gp$lwd  <- rep(gp$lwd,  length(polygons))
  if (length(gp$col) == 1)  gp$col  <- rep(gp$col,  length(polygons))
  if (length(gp$fill) == 1) gp$fill <- rep(gp$fill, length(polygons))
  
  for (i in seq_along(polygons)) {
    if (nrow(polygons[[i]]) >= 2) {
      grid::grid.polygon(
        polygons[[i]][, 1], polygons[[i]][, 2],
        default.units = "native",
        gp = grid::gpar(
          fill = gp$fill[i],
          col  = gp$col[i],
          lty  = gp$lty[i],
          lwd  = gp$lwd[i],
          lineend  = "butt",
          linejoin = "mitre"
        )
      )
    }
  }
  
  grid::seekViewport(name = paste0("hilbert_curve_", HilbertCurve:::get_plot_index(), "_global"))
  grid::upViewport()
  invisible(NULL)
}