#' Create Genomic Screenshots with BigWig Tracks and Gene Annotations
#'
#' @description
#' Creates publication-quality genomic screenshots displaying BigWig signal tracks,
#' BED regions, and optional gene annotations for specified genomic regions. Multiple
#' regions can be displayed side by side for comparison.
#'
#' @param bed Input genomic ranges, in any format compatible with ?importBed().
#' @param tracks Character vector of paths to visualization tracks:
#' \itemize{
#'   \item BigWig (.bw) files for continuous signal tracks
#'   \item BED files for discrete regions (any format compatible with importBed)
#' }
#' @param track.names Function or character vector for track labels. Default extracts names
#'   from file basenames.
#' @param col Colors for tracks. Default is a grey gradient from grey60 to grey10.
#' @param nbins The number of bins used to quantify the signal (overrides bw.n.breaks). Default= NA (no binning). 
#' @param bw.max Maximum value at which bigwig coverage values will be clipped. Default= NA (no clipping).
#' @param genome Genome assembly (e.g., "mm10", "dm6") for gene annotations.
#'   If specified, nearby genes will be displayed.
#' @param gtf Character path to a custom GTF file for custom gene annotations (alternative to the 'genome' parameter).
#' @param sel.gene.symbols If specified, only selected gene symbols will be plotted. Default= NULL.
#' @param border.col The color used for track borders. Default= NA.
#' @param border.lwd The line width used for track borders. Default= 1.
#' @param bw.min Minimum signal value at which bigwig coverage values will be clipped. Default= NA (no clipping).
#' @param bw.n.breaks If nbins is set to NA, specifies the maximum number of levels to which the signal
#' should be simplified to avoid polygons with too many points. Default= 100.
#' @param ngenes Integer. Number of nearest genes to display. Default = 1.
#' @param cex.gene.symbol Size of gene symbols. Default = 1.
#' @param offset.gene.symbol Vertical offset for gene symbols. Default= 1.
#' @param col.gene.ps Color for positive strand genes. Default = "tomato".
#' @param col.gene.ns Color for negative strand genes. Default = "cornflowerblue".
#' @param exon.border Color for exon rectangle borders. Default = NA, meaning no border.
#' @param gtf.transcript If specified, only features with the corresponding biotype will be plotted as transcripts. Default = NULL.
#' @param gtf.exon Name of exon features in the GTF file (default = "exon").
#' @param gtf.symbol Name of the attribute containing gene symbols in the GTF file. Default = "gene_symbol".
#' @param gtf.transcript.id Name of the attribute containing transcript IDs in the GTF file. Default = "transcript_id".
#' @param region.width Width of each genomic region. Default= 1.
#' @param space.width Horizontal space between regions. Default= 1.
#' @param bw.height Height of BigWig tracks. Default= 1.
#' @param bed.height Height of BED tracks. Default= 1.
#' @param space.height Vertical space between tracks. Default= 1.
#' @param gene.height Height of gene tracks. Default= 1.
#' @param gene.space.height Space between genes. Default= 1.
#' @param add Should the plot beadded to the existing one? Default= FALSE.
#'
#' @return
#' A screenshot of the genomic tracks/regions.
#'
#' @examples
#' # Plotting parameters for a nice layout
#' par(mai = c(0.9, 0.9, 0.9, 0.9),
#'     las = 1,
#'     tcl = -0.1,
#'     mgp = c(1.5, 0.35, 0),
#'     cex = 1,
#'     cex.lab = 9/12,
#'     cex.axis = 7/12,
#'     bty = "n",
#'     lend = 2)
#'
#' # Simple example with two regions and 4 tracks
#' bwScreenshot(
#'   bed= c("chr3R:30760926-30794202", "chr3R:30765926-30789202"),
#'   tracks= c("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PH18_merge.bw",
#'             "/groups/stark/vloubiere/projects/epigenetic_cancer/db/peaks/ATAC/ATAC_PH18_conf_peaks.narrowPeak",
#'             "/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PHD11_merge.bw",
#'             "/groups/stark/vloubiere/projects/epigenetic_cancer/db/peaks/ATAC/ATAC_PHD11_conf_peaks.narrowPeak"),
#'   genome= "dm6"
#' )
#'
#' # Genes only
#' bwScreenshot(
#'   bed= c("chr3R:30760926-30794202", "chr3R:30765926-30789202"),
#'   genome= "dm6"
#' )
#'
#' @return A genomic screenshot.
#' @export
bwScreenshot <- function(
    bed,
    tracks= character(),
    track.names= function(x) gsub("(.*)[.].*$", "\\1", x),
    col= colorRampPalette(c("grey60", "grey10"))(length(tracks)),
    bw.max= NA,
    nbins= 500,
    genome,
    gtf,
    sel.gene.symbols= NULL,
    border.col= NA,
    border.lwd= 1,
    bw.min= NA,
    bw.n.breaks= 100,
    ngenes= 1,
    cex.gene.symbol= 1,
    offset.gene.symbol= 0.25,
    col.gene.ps= "tomato",
    col.gene.ns= "cornflowerblue",
    exon.border= NA,
    gtf.transcript= NULL,
    gtf.exon= "exon",
    gtf.symbol= "gene_symbol",
    gtf.transcript.id= "transcript_id",
    region.width= 1,
    space.width= 1,
    bw.height= 1,
    bed.height= 1,
    space.height= 1,
    gene.height= 1,
    gene.space.height= 1,
    add= FALSE
)
{
  # Scale expansion factors ----
  border.lwd <- border.lwd*.5
  cex.gene.symbol <- cex.gene.symbol*.7
  region.width <- region.width*100
  space.width <- space.width*30
  bw.height <- bw.height*100
  bed.height <- bed.height*20
  space.height <- space.height*2
  gene.height <- gene.height*4
  gene.space.height <- gene.space.height*7
  
  # Import bed regions ----
  regions <- vlite::importBed(bed = bed)[, .(seqnames, start, end)]
  stopifnot(is.na(nbins) || min(regions[, end-start+1])>nbins)
  # Add index, width and plot limits
  regions[, region.idx:= .I]
  regions[, width:= end-start+1]
  regions[, xleft:= (region.idx-1)*(region.width+space.width)]
  regions[, xright:= xleft+region.width]
  
  # Import gene features ----
  if((!missing(genome) | !missing(gtf)) && ngenes>0)
  {
    genes <- .genomeGTFfeatures(
      genome= genome,
      regions= regions,
      ngenes= ngenes,
      sel.gene.symbols= sel.gene.symbols,
      gene.height= gene.height,
      gene.space.height= gene.space.height,
      gtf= gtf,
      gtf.transcript= gtf.transcript,
      gtf.exon= gtf.exon,
      gtf.symbol= gtf.symbol,
      gtf.transcript.id= gtf.transcript.id
    )
  }
  
  # Format tracks metadata table ----
  meta <- data.table(track.file= tracks,
                     track.class= .checkTrackClass(tracks),
                     track.col= col)
  if(nrow(meta))
  {
    # Add names
    if(is.function(track.names))
      track.names <- track.names(basename(tracks))
    meta[, track.name:= track.names]
    # Add index, track height and score cutoffs
    meta[, track.idx:= .I]
    meta[, track.height:= switch(track.class,
                                 "bw"= bw.height,
                                 "bed"= bed.height), track.class]
    meta[track.class=="bw", c("track.cutoff.min", "track.cutoff.max"):= .(bw.min, bw.max)]
    # Add plot limits
    meta[, ytop:= rev(cumsum(rev(track.height+space.height)))]
    if(exists("genes", inherits= FALSE) && nrow(genes)>0) # Adjust if genes are plotted
      meta[, ytop:= ytop+max(genes$ytop)+gene.space.height]
    meta[, ybottom:= ytop-track.height]
  }
  
  # Initiate plot ----
  if(!add)
  {
    # y maximum limit
    y.max <- if(nrow(meta))
      max(meta$ytop) else if(exists("genes", inherits= F))
        max(genes$ytop) else
          0
    # Empty plot
    plot(x= range(regions[, c(xleft, xright)]),
         y= c(0, y.max),
         type= "n",
         frame= FALSE,
         axes= FALSE,
         ylab= NA,
         xlab= "Genomic coordinates")
  }
  
  # Add regions coordinates on x axis ----
  regions[, {
    .name <- paste0(seqnames, ":",
                    formatC(start, big.mark = ","), "-",
                    formatC(end, big.mark = ","))
    axis(1,
         at = (xleft+xright)/2,
         labels = .name,
         tick = FALSE)
  }, (regions)]
  
  # Plot tracks if specified ----
  if(nrow(meta))
    meta[, {
      if(track.class=="bw")
      {
        # Plot method for bw
        .screenshotBwMethod(
          regions= regions,
          track.file= track.file[1],
          nbins= nbins,
          bw.n.breaks= bw.n.breaks,
          track.name= track.name[1],
          track.col= track.col[1],
          track.cutoff.min= track.cutoff.min[1],
          track.cutoff.max= track.cutoff.max[1],
          track.height= track.height[1],
          ybottom= ybottom[1],
          ytop= ytop[1],
          border.col= border.col,
          border.lwd= border.lwd
        )
      }else if(track.class=="bed")
      {
        # Plot method for bed
        .screenshotBedMethod(
          regions= regions,
          track.file= track.file[1],
          track.col= track.col[1],
          track.height= track.height[1],
          track.name= track.name[1],
          ybottom= ybottom[1],
          ytop= ytop[1],
          border.col= border.col,
          border.lwd= border.lwd
        )
      }
      .SD
    }, track.idx]
  
  # Plot scale bar ----
  regions[, {
    # Compute ideal scale bar
    bar <- 10^floor(log10(width))
    x0 <- xright-(bar/width*(xright-xleft))
    segments(x0,
             par("usr")[4],
             xright[1],
             par("usr")[4],
             xpd= T)
    # Simplif label
    bar <- if(bar>1e3)
      paste0(bar/1000, "kb") else if(bar>1e6)
        paste0(bar/1000, "Mb") else
          paste(bar, "bp")
    text((x0+xright)/2,
         par("usr")[4],
         bar,
         pos= 3,
         xpd= T)
  }, .(width, xright)]
  
  
  # Plot genes if genome specified ----
  if(exists("genes", inherits = F) && nrow(genes)>0)
  {
    .screenshotGtfMethod(
      regions= regions,
      genes= genes,
      col.ps= col.gene.ps,
      col.ns= col.gene.ns,
      exon.border= exon.border,
      offset.symbol= offset.gene.symbol,
      cex.symbol= cex.gene.symbol
    )
  }
}
