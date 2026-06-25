# Check the class (bw or bed) of the tracks provided int the ?screenshotBw() function
.checkTrackClass <- function(file)
{
  if(!is.character(file))
    stop("file should be character vector")
  type <- fcase(grepl(".bed$|.narrowPeak$|.broadPeak$", file), "bed",
                grepl(".bw$|.bigwig$|.bigWig$", file), "bw",
                grepl(".gtf$", file), "gtf",
                default = as.character(NA))
  if(anyNA(type))
    stop("File extension should be one of .bed, .narrowPeak, .broadPeak, .bw, .bigWig, .gtf")
  
  return(type)
}

# A set of parameters for plotting specific genomes using ?bwScreenshot()
.genomeGTFparameters <- function(genome)
{
  # gtf= valid path to a gtf file
  # gtf.transcript= transcript type name in the gtf file
  # gtf.exon= exon type name in the gtf file
  # gtf.symbol= column containing gene symbols in the gtf file
  # gtf.transcript.id= column containing transcript id in the gtf file
  params <- if(genome=="dm6")
  {
    list(gtf= "/zssd/scratch/vincent.loubiere/genomes/Drosophila_melanogaster/flybase/dm6/dmel-all-r6.36.gtf",
         gtf.transcript= NULL,
         gtf.exon= "exon",
         gtf.symbol= "gene_symbol",
         gtf.transcript.id= "transcript_id")
  } else if(genome=="mm10") {
    list(gtf= "/zssd/scratch/vincent.loubiere/genomes/Mus_musculus/GENCODE/gencode.vM25.basic.annotation_for_viewer.gtf.gz",
         gtf.transcript= NULL,
         gtf.exon= "exon",
         gtf.symbol= "gene_name",
         gtf.transcript.id= "transcript_id")
  } else
    stop("Genome not supported")
  
  # Return parameters
  return(params)
}

# Method to import gene features from gtf for plotting with ?bwScreenshot()
.genomeGTFfeatures <- function(
    genome,
    regions,
    ngenes,
    sel.gene.symbols= NULL,
    gene.height,
    gene.space.height,
    gtf,
    gtf.transcript,
    gtf.exon,
    gtf.symbol,
    gtf.transcript.id
)
{
  # Import parameters using helper function
  if(!missing(genome))
  {
    params <- .genomeGTFparameters(genome)
    list2env(params, envir = environment())
  }
  
  # Import features
  feat <- rtracklayer::import(
    gtf,
    colnames = c("type", gtf.symbol, gtf.transcript.id)
  )
  GenomeInfoDb::seqlevelsStyle(feat) <- "UCSC"
  feat <- vlite::importBed(feat)
  
  # Make column names and types reproducible
  setnames(feat, c(gtf.symbol, gtf.transcript.id), c("symbol", "id"))
  feat[, type:= as.character(type)]
  
  # Check
  if(!gtf.exon %in% feat$type)
    warning("gtf.exon type is missing from the provided gtf file. No genes will be plotted.")
  
  # Remove transcripts with no id
  if(anyNA(feat$id)) {
    feat[, {
      print(
        paste(
          sum(is.na(id)), "/", .N, shQuote(type), "gtf entries with NA", gtf.transcript.id, "were removed."
        )
      )
    }, type]
    feat <- feat[!is.na(id)]
    if(length(unique(feat$type))<2)
      warning("GTF: only one feature type was left after removing entries with NA.")
  }
  
  # Select specific transcripts if relevant
  if(!is.null(gtf.transcript)) {
    if(!gtf.transcript %in% feat$type)
      warning("gtf.transcript type is missing from the provided gtf file. No genes will be plotted.")
    feat <- feat[type %in% c(gtf.transcript, gtf.exon)]
    # Remove exons of non-selected transcripts
    feat <- feat[id %in% feat[type!=gtf.exon, id]]
  }
  
  # Clip genes to regions of interest
  feat <- vlite::clipBed(feat, regions)
  
  # Simplify
  feat[, type:= ifelse(type==gtf.exon, "exon", "transcript")]
  feat[type=="transcript", c("start", "end"):= .(min(start), max(end)), id]
  feat <- unique(feat)
  
  # Compute y order for plotting
  feat[, ov:= collapseBed(feat, return.idx.only = TRUE)]
  setorderv(feat, "start")
  feat[, y:= as.numeric(factor(id, unique(id))), ov]
  
  # Select n first genes to be plotted
  feat <- feat[y <= ngenes]
  
  # Select gene symbols of interest
  if(!is.null(sel.gene.symbols))
    feat[!symbol %in% sel.gene.symbols, symbol:= NA]
  
  # Add plot limits
  feat[, ytop:= (max(c(y, 0))-y+1)*(gene.height+gene.space.height)]
  feat[, ybottom:= ytop-gene.height]
  
  # Return gene features (to plot with .screenshotGtfMethod)
  return(feat)
}

# Method to plot bw file in bwScreenshot ----
.screenshotBwMethod <- function(
    regions,
    track.file,
    nbins,
    track.name,
    track.col,
    track.cutoff.min,
    track.cutoff.max,
    show.bw.range,
    track.height,
    ybottom,
    ytop,
    border.col,
    border.lwd
)
{
  if(is.null(nbins)) {
    # Import bw score
    gr <- GenomicRanges::GRanges(regions)
    sel <- rtracklayer::BigWigSelection(gr, "score")
    var <- rtracklayer::import.bw(track.file, selection= sel)
    var <- data.table::as.data.table(var)
    var <- var[, .(seqnames, start, end, score)]
    # Replace NAs with 0
    var[is.na(score), score:= 0]
    # Fill gaps with 0s
    gaps <- subtractBed(regions[, .(seqnames, start, end)],
                        var)
    gaps[, score:= as.numeric(0)]
    if(nrow(gaps))
      var <- rbind(var, gaps, fill = TRUE)
    # Extract bins for each region
    var <- regions[, {
      # Clip Polygons to region
      clipBed(var, .SD)
    }, .(region.idx, width, xleft, xright)]
  } else {
    # Bin regions
    var <- binBed(regions, nbins = nbins)
    var[, score:= bwCoverage(var, track.file)]
    var[is.na(score), score:= 0]
    var <- var[, .(region.idx, width, xleft, xright, seqnames, start, end, score)]
  }
  # Clip score based on cutoffs
  if(is.na(track.cutoff.min))
    track.cutoff.min <- min(c(0, var$score))
  if(is.na(track.cutoff.max))
    track.cutoff.max <- max(var$score)
  var[score<track.cutoff.min, score:= track.cutoff.min]
  var[score>track.cutoff.max, score:= track.cutoff.max]
  # Compute range
  score.range <- track.cutoff.max-track.cutoff.min
  # Scale signal for plotting
  var[, ypos:= (score-track.cutoff.min)/score.range]  # From 0 to 1
  var[, ypos:= ypos*track.height] # Scale to track height
  var[, ypos:= ypos+ybottom] # align to bottom
  # Baseline
  baseline <- -track.cutoff.min/score.range*track.height+ybottom
  # Plot track name (y axis)
  axis(2,
       at= ybottom+track.height/2,
       labels = track.name,
       lwd= 0,
       las= 1,
       cex.axis= par("cex.lab"))
  # Plot range limits (y axis)
  if(show.bw.range) {
    axis(2,
         at= ytop,
         labels = formatC(track.cutoff.max,
                          format = "e",
                          digits = 1),
         lwd= 0,
         las= 1)
    if(track.cutoff.min != 0) {
      axis(2,
           at= ybottom,
           labels = formatC(track.cutoff.min,
                            format = "e",
                            digits = 1),
           lwd= 0,
           las= 1)
    }
  }
  # Scale polygons and plot ----
  var[, coor:= rowMeans(.SD), .SDcols= c("start", "end")]
  var[, coor:= (coor-start[1])/width*diff(c(xleft, xright))+xleft, .(region.idx, width, xleft, xright)]
  # plot
  var[, {
    polygon(
      c(xleft, coor, xright),
      c(baseline, ypos, baseline),
      col= track.col,
      border= border.col,
      lwd= border.lwd
    )
  }, .(region.idx, xleft, xright)]
}

# Method to plot bed file in bwScreenshot----
.screenshotBedMethod <- function(
    regions,
    track.file,
    track.col,
    track.height,
    track.name,
    ybottom,
    ytop,
    border.col,
    border.lwd
)
{
  # Import bed
  var <- importBed(track.file)
  var[, score:= 1]
  # Check with of regions to plot
  if(any(var[,end-start]==0) & is.na(border.col))
    warning("Some regions in bed file(s) are width 1 and might not appear because border.col is set to NA. Consider border.col= 'black'")
  # Plot
  regions[, {
    # Clip Polygons to region
    poly <- clipBed(var, .SD[, .(seqnames, start, end)])
    if(nrow(poly)) {
      # Adjust start and end
      poly$start <- (poly$start-start)/width*diff(c(xleft, xright))+xleft
      poly$end <- (poly$end-start)/width*diff(c(xleft, xright))+xleft
      # Plot polygons
      poly[, {
        rect(start,
             ybottom,
             end,
             ytop,
             col= track.col,
             border= border.col,
             lwd = border.lwd)
      }]
    }
  }, region.idx]
  # Plot label (y axis)
  axis(2,
       at= ybottom+track.height/2,
       labels = track.name,
       lwd= 0,
       las= 1,
       cex.axis= par("cex.lab"))
}

# Method to plot gene from GTF in genesScreenshot ----
.screenshotGtfMethod <- function(regions,
                                 genes,
                                 col.ps,
                                 col.ns,
                                 exon.border,
                                 offset.symbol,
                                 cex.symbol)
{
  regions[, {
    # Select transcripts overlapping current region
    current <- clipBed(genes, .SD)
    if(nrow(current)) {
      # Scale start and end
      current$start <- (current$start-start)/width*diff(c(xleft, xright))+xleft
      current$end <- (current$end-start)/width*diff(c(xleft, xright))+xleft
      # Only keep the gene symbol for the lowest transcript on the plot
      current[type=="transcript", symbol:= ifelse(y<max(y), NA, symbol), symbol]
      # Plot gene bodies (segments)
      current[type=="transcript", {
        segments(start,
                 (ytop+ybottom)/2,
                 end,
                 (ytop+ybottom)/2,
                 col= ifelse(strand=="-", col.ns, col.ps))
        if(.N>0)
          text((start+end)/2,
               ybottom,
               symbol,
               pos= 1,
               offset= offset.symbol,
               cex= cex.symbol,
               xpd= NA)
      }]
      # Plot exons (rectangles)
      current[type=="exon", {
        rect(start,
             ybottom,
             end,
             ytop,
             col = ifelse(strand=="-", col.ns, col.ps),
             border= exon.border)
      }]
    }
  }, region.idx]
}
