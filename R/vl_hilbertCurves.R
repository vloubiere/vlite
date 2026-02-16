#' Function to plot hilbert curves
#' 
#' Wrapper around ?HilbertCurve package (mode= "pixel").
#'
#' @param track Path to the bw file to be plotted.
#' @param title The tile of the heatmap. Default= "pixel mode".
#' @param genome The genome to bin. Default= "dm6.
#' @param width The width of the bins used to tile the genome. Default= 2500.
#' @param canonical.chr The canonical chromosomes that should be plotted. Default= c("chrX", "chr2L", "chr2R", "chr3L", "chr3R", "chr4").
#' @param clip The values used to clip the data. Default= quantile(x, .9995).
#' @param highlight A set of regions that should be higlighted on the heatmap. The content of the 'name' column will be added as text.
#'
#' @returns
#' @export
#'
#' @examples
vl_hilbertCurves <- function(
    track,
    title= "pixel mode",
    genome= "dm6",
    width= 2500,
    canonical.chr= c("chrX", "chr2L", "chr2R", "chr3L", "chr3R", "chr4"),
    clip= quantile(cov, .9995),
    highlight= NULL
)
{
  # Bin genome ----
  chr <- vlite::getBSgenomeSize(genome = genome)
  chr <- chr[seqnames %in% canonical.chr]
  bins <- binBed(
    chr,
    bins.width = width
  )
  bins[, idx:= .I]
  # Compute coverage and clip extremes ----
  cov <- vlite::bwCoverage(bins, track)
  cov[cov>clip] <- clip
  print(paste("Max value set to", clip))
  # Color function ----
  col_fun = circlize::colorRamp2(c(0, max(cov)), c("white", "red"))
  # Initiate curve ----
  hc = HilbertCurve::HilbertCurve(
    1,
    max(bins$idx)+1,
    level = 10,
    mode = "pixel",
    title = title,
    legend= FALSE
  )
  # Add signal layer ----
  x1 <- bins$idx
  x2 <- bins$idx+1
  HilbertCurve::hc_layer(
    object = hc,
    x1 = x1,
    x2 = x2,
    mean_mode = "absolute",
    col = col_fun(cov)
  )
  # Add chromosome polygons ----
  x1 <- bins[, idx[1], seqnames]$V1
  x2 <- bins[, idx[.N], seqnames]$V1+1
  HilbertCurve::hc_text(
    object = hc,
    x1 = x1,
    x2 = x2, 
    labels = unique(bins$seqnames),
    gp = grid::gpar(fontsize = 10),
    centered_by = "polygon"
  )
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # The current HilbertCurve::hc_polygon function is bugged (legend= FALSE is there twice)
  # So for now I overwrote it with hc_polygon_fixed
  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  hc_polygon_fixed(
    object = hc,
    x1 = x1,
    x2 = x2
  )
  # highlight genes ----
  if(!is.null(highlight))
  {
    highlight <- importBed(highlight)
    highlight[, {
      # Polygons ----
      .c <- intersectBed(bins, highlight)
      hc_polygon_fixed(
        object = hc,
        x1 = .c[1,idx],
        x2 = .c[.N,idx]+1
      )
      HilbertCurve::hc_text(
        object = hc,
        x1 = .c[1,idx],
        x2 = .c[.N,idx],
        labels = paste0("   ", unique(name)),
        gp = grid::gpar(fontsize = 10),
        centered_by = "polygon",
        just= "left"
      )
    }, (highlight)]
  }
}