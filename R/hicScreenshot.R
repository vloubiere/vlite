#' Function to plot a HiC map with bw tracks on the bottom
#' 
#' Wrapper around HiCExperiment::import and vlite::bwScreenshot to plot a HiC matrix with bw tracks.
#'
#' @param mcool.file Path to the mcool file to be plotted.
#' @param mcool.file Path to the misha file to be plotted.
#' @param region String specifying a unique region to plot, in the format "chr2L:16300000-16600000".
#' @param resolution The resolution to be used. Only the ones existing in the mcool file will work (see ?HiCExperiment::availableResolutions)
#' @param map.name The name of the hic map. By default, the basename of the mcool file will be used.
#' @param pdf.file Output pdf file. Should automatically adapts to the tracks being plotted.
#' @param screenshot.cex.height Expansion factor to apply to the height of the bottom screenshot. Default= 1.
#' For all other arguments, please refer to ?vlite::bwScreenshot.
#'
#' @returns A plot of a HiC matrix with specific tracks
#' @export
#'
#' @examples
hicScreenshot <- function(
    mcool.file,
    misha.track,
    region,
    resolution,
    map.name= gsub(".mcool$|.track$", "", basename(mcool.file)),
    pdf.file= NULL,
    screenshot.cex.height= 1,
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
) {
  # Check region format ----
  region <- importBed(region)
  region.width <- region[, end-start+1]
  region <- region[, paste0(seqnames, ":", start, "-", end)]
  if(length(region)!=1)
    stop("For now hicScreenshot can only plot one region at once. Refer to you bioinformatics overlord bitteschönn")
  stopifnot(missing(mcool.file) | missing(misha.track))
  
  # Dispatch based on the input file ----
  if(!missing(mcool.file)) {
    
    # Import hic data ----
    hic <- HiCExperiment::import(
      con = mcool.file,
      resolution = resolution,
      focus = region
    )
    dmat <- as.data.table(hic@interactions)
    dmat[, score:= hic@scores$balanced]
  } else if(!missing(misha.track)) {
    
    # Load misha db ----
    library("misha", lib.loc = "/home/michael.szalay/anaconda3/envs/misha/lib/R/library/")
    # library("misha")
    mDBloc <- '/zdata/data/mishaDB/trackdb/'
    db <- 'dm6'
    dbDir <- paste0(mDBloc, db, '/')
    gdb.init(dbDir)
    gdb.reload()
    misha.track <- gtrack.ls(misha.track)
    
    # Import misha track data ----
    coor <- unlist(tstrsplit(region, ":|-"))
    chrom <- coor[1]
    start <- as.integer(coor[2])
    end <- as.integer(coor[3])
    
    # Create iterator ----
    interval2D <- gintervals.2d(chrom, start, end, chrom, start, end)
    binnedIterator <- giterator.intervals(
      intervals= interval2D,
      iterator= c(resolution, resolution)
    )
    
    # Extract iterator data ----
    data <- gextract(
      misha.track,
      binnedIterator,
      iterator= binnedIterator,
      colnames = "score"
    )
    
    # Format similar to cool files ----
    dmat <- as.data.table(data)
    colnames(dmat) <- c(
      "seqnames1", "start1", "end1",
      "seqnames2", "start2", "end2",
      "score", "intervalID"
    )
    dmat$bin_id1 <- (dmat$end1-start) / resolution
    dmat$bin_id2 <- (dmat$end2-start) / resolution
  } else 
    stop("input file could not be determined.")
  
  # Adjust region to closest hic bins ----
  adj.region <- paste0(dmat[1, seqnames1], ":", dmat[1, start1], "-", dmat[.N, end2])
  if(adj.region!=region)
    message(paste0("Target region was adjusted to ", adj.region))
  
  # Fill gaps ----
  bin.idx.1 <- range(dmat[, bin_id1])
  bin.idx.2 <- range(dmat[, bin_id2])
  dmat <- merge(
    CJ(bin_id1= bin.idx.1[1]:bin.idx.1[2], bin_id2= bin.idx.2[1]:bin.idx.2[2]),
    dmat[, .(bin_id1, bin_id2, score)],
    all.x= T
  )
  
  # HiC matrix ----
  mat <- dcast(dmat, bin_id1~bin_id2, value.var = "score")
  mat <- as.matrix(mat, 1)
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  
  # Set plotting parameters ----
  if(!missing(mcool.file)) {
    # Color
    Cc <- c(
      "#FFFFFF",
      "#FFFFCC",
      "#FFEDA0",
      "#FED976",
      "#FEB24C",
      "#FD8D3C",
      "#FC4E2A",
      "#E31A1C",
      "#BD0026",
      "#800026",
      "#000000"
    )
    Cc <- colorRampPalette(Cc)(256)
    
    # Min and max values are defined on a log scale
    vmin <- 0.0001
    vmax <- 0.1
    mat_log <- log10(mat)
    breaks <- seq(
      log10(vmin),
      log10(vmax),
      length.out = length(Cc) + 1
    )
    
  } else {
    
    # Color
    Cc <- c("darkblue", "white", "darkred")
    Cc <- colorRampPalette(Cc)(200)
    
    # Min and max values 
    vmin <- -100
    vmax <- 100
    breaks <- seq(
      vmin,
      vmax,
      length.out= length(Cc) + 1
    )
  }
  
  # Initiate plot ----
  if(is.character(pdf.file)) {
    # Top and bottom margins
    top.margin.in <- .75
    bottom.margin.in <- .5
    # bottom screenshot height
    bw.height.cm <- bw.height*1.5*sum(grepl(".bw$", tracks))
    bed.height.cm <- bed.height*0.3*sum(grepl(".bed$", tracks))
    space.height.cm <- space.height*0.03*(length(tracks)-1)
    gene.height.cm <- gene.height*0.06*(!missing(genome) | !missing(gtf))*ngenes
    gene.space.height.cm <- gene.space.height*0.105*(ngenes-1)
    screenshot.height.cm <- (bw.height.cm+bed.height.cm+space.height.cm+gene.height.cm+gene.space.height.cm)*(screenshot.cex.height*.7)
    # Open pdf
    pdf(
      pdf.file,
      width  = 20 / 2.54,
      height = (12 + screenshot.height.cm) / 2.54 + top.margin.in + bottom.margin.in
    )
    # Fixed layout
    layout(
      matrix(c(1, 2), nrow = 2),
      widths  = lcm(12),
      heights = c(lcm(12), lcm(screenshot.height.cm)),
      respect = TRUE
    )
    vl_par(
      mai = c(0, 0, 0, 0),
      omi = c(bottom.margin.in, 0, top.margin.in, 0),
      xaxs = "i",
      yaxs = "i",
      xpd= NA
    )
  }
  
  # Plot heatmap ----
  vl_heatmap(
    if(!missing(mcool.file)) mat_log else mat,
    breaks= breaks,
    cluster.rows = F,
    show.rownames = F,
    show.colnames = F,
    col= Cc,
    na.col = "white",
    show.legend = !missing(misha.track),
    legend.title = "Score"
  )
  title(
    main= map.name,
    outer = T,
    line = 2
  )
  
  # Add special heatkey ----
  if(!missing(mcool.file)) {
    vlite::heatkey(
      breaks= seq(log10(vmin), log10(vmax), length.out= 256),
      labels = seq(log10(vmin), log10(vmax), 1),
      log10.labels= T, # 10^labels will be plotted instead of actual labels
      col= Cc, 
      main = "ICE norm."
    )
    # Add scale bar
    bar <- 10^floor(log10(region.width))
    xleft <- par("usr")[1]
    xright <- par("usr")[2]
    x0 <- xright-(bar/region.width*(xright-xleft))
    scale.y <- par("usr")[4]+strheight("M")/2
    segments(x0,
             scale.y,
             xright,
             scale.y,
             xpd= NA)
    # Simplif label
    bar <- if(bar>1e3)
      paste0(bar/1000, "kb") else if(bar>1e6)
        paste0(bar/1000, "Mb") else
          paste(bar, "bp")
    text((x0+xright)/2,
         scale.y,
         bar,
         pos= 3,
         offset= 0.15,
         xpd= NA)
  }
  
  # Plot tracks ----
  vlite::bwScreenshot(
    bed= adj.region,
    tracks= tracks,
    track.names= track.names,
    col= col,
    bw.max= bw.max,
    nbins= nbins,
    genome= genome,
    gtf= gtf,
    sel.gene.symbols= sel.gene.symbols,
    border.col= border.col,
    border.lwd= border.lwd,
    bw.min= bw.min,
    ngenes= ngenes,
    cex.gene.symbol= cex.gene.symbol,
    offset.gene.symbol= offset.gene.symbol,
    col.gene.ps= col.gene.ps,
    col.gene.ns= col.gene.ns,
    exon.border= exon.border,
    gtf.transcript= gtf.transcript,
    gtf.exon= gtf.exon,
    gtf.symbol= gtf.symbol,
    gtf.transcript.id= gtf.transcript.id,
    region.width= region.width,
    space.width= space.width,
    bw.height= bw.height,
    bed.height= bed.height,
    space.height= space.height,
    gene.height= gene.height,
    gene.space.height= gene.space.height,
    show.scale.bar= F,
    add= add
  )
  
  if(is.character(pdf.file))
    dev.off()
}