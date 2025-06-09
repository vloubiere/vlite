#' Import Genomic Ranges in BED Format
#'
#' @description
#' Imports genomic coordinates from various input formats into a standardized data.table format.
#' Supported input types: GRanges objects, character strings ("chr2L:30000-31000:-"), .bed, .narrowPeak
#' or .broadPeak file paths (see details).
#'
#' @param bed Input genomic ranges in one of the following formats:
#' \itemize{
#'   \item Character vector of genomic coordinates ("chr:start-end[:strand]").
#'   \item File path to .bed, .narrowPeak, or .broadPeak file.
#'   \item GRanges object.
#'   \item data.frame or data.table with required columns.
#' }
#' @param col.names Can be used to overwrite default column names.
#'
#' @details
#' Character coordinates should follow the syntax "chr:start-end[:strand]".
#' File specifications can be found at https://genome.ucsc.edu/FAQ/FAQformat.html.
#' GRanges objects are directly converted to data.table and metadata columns are preserved.
#' Required columns for data.frames/data.tables are 'seqnames', 'start', 'end' and 'strand'.
#' If missing, end and strand will be set to 'start' and '*', respectively.
#'
#' @examples
#' # Import from character coordinates
#' importBed(c("chr2L:1000-2000:+", "chr3L:1000-2000:-", "chr3L:4000-6000:*"))[]
#'
#' # From GRanges
#' importBed(GRanges("chr2L", IRanges(c(1000, 5000), c(2000, 6000))))[]
#'
#' # From bed file
#' bed_file <- tempfile(fileext = ".bed")
#' exportBed(gr, bed_file)
#' importBed(bed_file)[]
#'
#' # From narrowpeak file
#' peaks <- data.table(
#'   seqnames = "chr2L",
#'   start = 1000,
#'   end = 2000,
#'   name = "peak1",
#'   score = 100,
#'   strand = "+",
#'   signalValue = 5.5,
#'   pValue = 0.001,
#'   qValue = 0.05,
#'   peak = 1500
#' )
#' narrowPeak_file <- tempfile(fileext = ".narrowPeak")
#' importBed(narrowPeak_file)
#'
#' @export
importBed <- function(bed,
                      col.names= NULL)
{
  if(is.character(bed)) {
    if(any(grepl("\\.(bed|narrowPeak|broadPeak)$", bed))) {

      # Import bed file using rtracklayer ----
      current <- as.data.table(rtracklayer::import(bed))

    } else {

      # Import character coordinates using GenomicRanges ----
      gr <- gsub(",", "", bed) # Remove potential comas
      gr <- GenomicRanges::GRanges(gr)
      check.col <- names(mcols(gr))
      current <- as.data.table(gr)
      # If no 'width' column in mcols, remove it
      if(!length(check.col) || !"width" %in% check.col)
        current$width <- NULL
      # Keep the provided coordinates as names
      current$name <- bed
    }
  } else if(is.data.frame(bed)) {

    # Coerce to data.table ----
    current <- if(is.data.table(bed)) data.table::copy(bed) else as.data.table(bed)
    # Check required columns
    if(!all(c("seqnames", "start") %in% names(current)))
      stop("Genomic ranges should contain at least 'seqnames' and 'start' fields.")
    # Format required columns
    current[, seqnames := as.character(seqnames)]
    current[, start := as.integer(start)]
    if(!"end" %in% names(current)) current[, end:= as.integer(start)]
    current[, end:= as.integer(end)]
    if(!"strand" %in% names(current)) current[, strand:= "*"]
    # Check strand values
    if(any(!current$strand %in% c("+", "-", "*")))
      stop("All strand values should be one of c('+', '-' or '*') -> malformed genomic ranges!")

  } else if(class(bed)[1]=="GRanges") {

    # Coerce to data.table ----
    gr <- GenomicRanges::GRanges(bed)
    check.col <- names(mcols(gr))
    current <- as.data.table(gr)
    # If no 'width' column in mcols, remove it
    if(!length(check.col) || !"width" %in% check.col)
      current$width <- NULL

  } else {
    stop("Input format could not be determined. See ?vlite::importBed.")
  }

  # Return bed data.table ----
  return(current)
}

#' Export Genomic Ranges to BED Format Files
#'
#' @description
#' Exports genomic ranges as .bed, .narrowPeak or .broadPeak files. See details.
#'
#' @param bed Input genomic ranges, in any format compatible with ?importBed.
#' @param file Output file path ending with .bed, .narrowPeak or .broadPeak extension.
#'
#' @details
#' File specifications can be found at https://genome.ucsc.edu/FAQ/FAQformat.html.
#' Upon saving, 1-based start coordinates will be converted to 0-based, following BED format specifications.
#'
#' @examples
#' # Export simple BED format:
#' exportBed("chr2L:1000-2000:+", file = "test.bed")
#' # Lookup saved file
#' fread("test.bed")
#'
#' # Export narrowPeak format:
#' peak_data <- data.table(
#'   seqnames = "chr2L",
#'   start = 1000,
#'   end = 2000,
#'   name = "peak1",
#'   score = 100,
#'   strand = "+",
#'   signalValue = 5.5,
#'   pValue = 0.001,
#'   qValue = 0.05,
#'   peak = 1500
#' )
#' exportBed(bed= peak_data, file = "test.narrowPeak")
#' # Lookup saved file
#' fread("test.narrowPeak")
#'
#' @export
exportBed <- function(bed, file)
{
  # Output type ----
  type <- gsub(".*[.](.*)$", "\\1", file)
  if(!type %in% c("bed", "narrowPeak", "broadPeak"))
    stop("Path should have .bed, .narrowPeak or .broadPeak extension.")

  if(type=="bed") {

    # Export using rtracklayer ----
    rtracklayer::export(GenomicRanges::GRanges(bed), file)

  } else {

    # Custom method for broadPeak/narrowPeak ----
    current <- vlite::importBed(bed)
    # Check if required columns exist
    col.names <- if(type=="broadPeak") {
      c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
    } else if(type=="narrowPeak") {
      c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
    }
    missing.cols <- setdiff(col.names, names(current))
    if(length(missing.cols))
      stop(paste("Missing columns for" , type, "format:", paste0(missing.cols, collapse = ", ")))
    # Select and order columns
    current <- current[, (col.names), with= FALSE]
    # Convert start to 0 base and make sure coor are integers
    current[, start:= start-1]
    # Set scipen to avoid scientific notation
    options(scipen = 999)
    # Save
    fwrite(current,
           file = file,
           col.names = FALSE,
           row.names = FALSE,
           sep= "\t",
           quote= FALSE,
           na= NA)

  }
}

#' Resize Genomic Regions in a BED File
#'
#' @description
#' Resize genomic ranges.
#'
#' @param bed Input regions in any format compatible with ?importBed.
#' @param center Specifies the origin for resizing. Options are:
#' \itemize{
#'   \item `center`: Resize symmetrically around the midpoint of the region (default).
#'   \item `start`: Resize from the start of the region.
#'   \item `end`: Resize from the end of the region.
#'   \item `region`: Extend or contract the entire region's boundaries.
#' }
#' @param upstream Size of the extension upstream of the specified origin. Default= 500L.
#' @param downstream Size of the extension downstream of the specified origin. Default= 500L.
#' @param genome A BS genome name (e.g. "dm6", "mm10") used to clipped resized regions.
#' Default= NULL (no clipping).
#' @param ignore.strand If set to TRUE, resizing always proceeds from the leftmost coordinates,
#' irrespective of the strand. If set to FALSE (default), upstream and downstream resizing respect
#' the strand.
#'
#' @return
#' A data.table containing resized regions.
#'
#' @examples
#' # Import example BED file
#' bed <- importBed(c("chr2L:1000-2000:+", "chr2L:1000-2000:-"))
#'
#' # Resize using different centering options
#' resizeBed(bed, center = "start", upstream = 100, downstream = 200)[]
#' resizeBed(bed, center = "center", upstream = 100, downstream = 200)[]
#' resizeBed(bed, center = "end", upstream = 100, downstream = 200)[]
#' resizeBed(bed, center = "region", upstream = 100, downstream = 200)[]
#'
#' # Ignore strand during resizing
#' resizeBed(bed, center = "center", upstream = 100, downstream = 200, ignore.strand = TRUE)[]
#'
#' @export
resizeBed <- function(bed,
                      center= "center",
                      upstream= 500L,
                      downstream= 500L,
                      genome= NULL,
                      ignore.strand= FALSE)
{
  # Import bed for incapsulation ----
  current <- vlite::importBed(bed)

  # Get region anchor depending on strand ----
  newStart <- if(center=="start"){
    current[, ifelse(strand=="-" & !ignore.strand, end, start)]
  } else if(center=="end") {
    current[, ifelse(strand=="-" & !ignore.strand, start, end)]
  } else if(center=="center") {
    current[, floor(rowMeans(.SD)), .SDcols= c("start", "end")]
  } else if(center=="region") {
    current[, start]
  }

  # Compute width if relevant ----
  if(center=="region")
    width <- current[, end-start] # In this case, do not subtract 1!

  # Resize ----
  current$start <- newStart-ifelse(current$strand=="-" & !ignore.strand, downstream, upstream)
  current$end <- newStart+ifelse(current$strand=="-" & !ignore.strand, upstream, downstream)
  if(center=="region")
    current$end <- current$end+width

  # Clip regions outside of chromosomes ----
  if(!is.null(genome)) {
    sizes <- vlite::getBSgenomeSize(genome = genome)
    current <- vlite::clipBed(current, sizes)
  }

  # Return ----
  return(current)
}

#' Bin Genomic Regions into Fixed-Size or Sliding Windows
#'
#' @description
#' A wrapper around ?GenomicRanges::tile and ?GenomicRanges::slidingWindows, that creates fixed-size
#' or sliding window bins from genomic ranges.
#'
#' @param bed Input genomic ranges in any format compatible with ?importBed.
#' @param nbins Integer specifying the number of equal-sized bins in which each region is divided.
#' @param center.nbins If set to TRUE and nbins is specified, bins are centered around the regions' midpoints.
#' If set to FALSE (default), the first bin will start from the leftmost coordinate.
#' @param genome When center.bins is set to TRUE, a BS genome name (e.g. 'dm6', 'mm10') that will be used to clip
#' resized regions that might extend beyond chromosome size (see details).
#' @param bins.width Integer specifying the width of the bins in which each region is divided.
#' Used only when nbins is not specified.
#' @param steps.width Integer specifying the distance between starts of consecutive bins. Used in
#' combination with bins.width. Default= bins.width, resulting in non-overlapping bins.
#' Smaller values create overlapping bins.
#' @param ignore.strand Although bins will always start from the leftmost coordinates, this argument
#' controls whether bins' indices respect feature's orientation. Default= FALSE (i.e. downstream bins
#' receive higher indices).
#'
#' @details
#' **Centered Binning**:
#' Centering bins (center.nbins = TRUE) is useful for analyzing features relative to central points.
#' In this case, each region will first be extended symmetrically by region_width/(nbins-1)/2, and a
#' genome can be specified to avoid regions outside of chromosome sizes.
#'
#' @return A gr data.table with columns:
#' \itemize{
#'   \item line.idx: line index of the region in bed.
#'   \item bin.idx: bin index (unique for each line.idx).
#'   \item seqnames: chromosome or sequence name.
#'   \item start: bin start position.
#'   \item end: bin end position.
#'   \item Additional columns within bed input are preserved.
#' }
#'
#' @examples
#' # Equal division binning
#' regions <- importBed(c("chr2L:1-10:+", "chr2L:100-200:-"))
#'
#' # Split into 5 equal bins
#' binBed(regions, nbins = 5)
#'
#' # Non-overlapping 50bp bins
#' binBed(regions, bins.width = 50, steps.width = 50)
#'
#' # Overlapping 50bp bins
#' binBed(regions, bins.width = 50, steps.width = 25)
#'
#' # Compare start-anchored vs. centered binning
#' nc1 <- binBed(regions, nbins= 5)
#' c1 <- binBed(regions, nbins= 5, center.nbins = TRUE)
#' nc2 <- binBed(regions, bins.width = 25)
#'
#' plot(x= c(50, 250), y= c(4, 10), type= "n", xlab= "coordinates", ylab= NA, yaxt= "n")
#' abline(v= c(100, 200))
#' nc1[line.idx==2][, {text(x= mean(c(start[1], end[.N])), y= 9, pos= 3, "nbins=5"); rect(start, 8, end, 9)}]
#' c1[line.idx==2][, {text(x= mean(c(start[1], end[.N])), y= 7, pos= 3, "centered"); rect(start, 6, end, 7)}]
#' nc2[line.idx==2][, {text(x= mean(c(start[1], end[.N])), y= 5, pos= 3, "bins.width=25"); rect(start, 4, end, 5)}]
#'
#' @export
binBed <- function(bed,
                   nbins,
                   center.nbins= FALSE,
                   genome= NULL,
                   bins.width = NULL,
                   steps.width = bins.width,
                   ignore.strand= FALSE)
{
  # Check method ----
  method <- if(!missing(nbins)) {
    if(nbins %% 1 != 0 | nbins<1) {
      stop("nbins should be an integer >= 1!")
    } else
      "nbins"
  } else if(is.numeric(bins.width) & is.numeric(steps.width)) {
    if(any(c(bins.width, steps.width) %% 1 != 0) || any(c(bins.width, steps.width)<1)) {
      stop("bins.width and steps.width should be integers >= 1!")
    } else
      "slidingWindow"
  } else
    stop("nbins or bins.width should be specified!")
  if(center.nbins && method!="nbins")
    stop("The 'center.nbins' option only available when 'nbins' is provided.")
  if(any(c("line.idx", "bin.idx") %in% names(bed)))
    warning("'line.idx' or 'bin.idx' column(s) already exist in the input bed and will be overwritten.")

  # Import bed for incapsulation ----
  bed <- vlite::importBed(bed)
  bed[, line.idx:= .I]

  # Resize if bins have to be centered ----
  if(center.nbins) {
    # Compute ext size
    ext <- floor(bed[, (end-start+1)]/(nbins-1)/2)
    # Extend
    bed <- vlite::resizeBed(bed,
                            center = "region",
                            upstream = ext,
                            downstream = ext,
                            genome= genome,
                            ignore.strand = ignore.strand)
  }

  # Compute bins ----
  bins <- if(method=="nbins") {
    GenomicRanges::tile(GenomicRanges::GRanges(bed),
                        n = nbins)
  } else if(method=="slidingWindow") {
    GenomicRanges::slidingWindows(GenomicRanges::GRanges(bed),
                                  width= bins.width,
                                  step = steps.width)
  }

  # Expand list ----
  cols <- setdiff(names(bed),
                  c("seqnames", "start", "end", "strand", "line.idx"))
  res <- bed[, {
    c(as.data.table(..bins[[line.idx]]), .SD)
  }, line.idx, .SDcols= cols]

  # Bin indices ----
  res[, bin.idx:= rowid(line.idx)]
  if(!ignore.strand)
    res[strand=="-", bin.idx:= rev(bin.idx), line.idx]

  # Order columns and return ----
  setcolorder(res, c("line.idx", "bin.idx"))
  return(res)
}

#' Collapse or Merge Overlapping Genomic Ranges
#'
#' @description
#' A wrapper around ?GenomicRanges::reduce that merges overlapping genomic ranges.
#'
#' @param bed Input regions in any format compatible with ?importBed.
#' @param min.gapwidth Minimum distance between features to be merged. Default= 1L.
#' \itemize{
#'   \item Positive value: Merges features separated by ≤ min.gapwidth bases.
#'   \item Zero: touching yet non-overlapping features will not be merged.
#' }
#' @param return.idx.only If set to TRUE, returns the indices of merged regions
#'   instead of collapsing them. Default= FALSE.
#' @param ignore.strand If set to FALSE, only overlapping regions that are on the same strand are merged.
#' If set to TRUE (default), regions are merged irrespective of their strand.
#'
#' @return
#' If return.idx.only = FALSE: a gr data.table with columns:
#' \itemize{
#'   \item seqnames: chromosome or sequence name.
#'   \item start: start position of merged region.
#'   \item end: end position of merged region.
#'   \item strand: strand (*= unstranded).
#' }
#' If return.idx.only = TRUE: a run-length type id indicating overlapping regions belong
#' that can be merged.
#'
#' @examples
#' # Create example regions
#' bed <- importBed(c("chr2R:1000-2000:+", "chr2R:2001-2500:+", "chr2R:2300-3000:-", "chr2R:3500-4000:+"))
#'
#' # Merge overlapping or touching regions
#' collapseBed(bed)
#'
#' # Only merge overlapping regions
#' collapseBed(bed,  min.gapwidth= 0)
#'
#' # Only merge regions that are on the same strand
#' collapseBed(bed, ignore.strand = FALSE)
#'
#' # Get merge indices instead of merging
#' collapseBed(bed, ignore.strand = FALSE, return.idx.only = TRUE)
#'
#' # Merge regions within 500bp of each other
#' collapseBed(bed, min.gapwidth = 500)
#' collapseBed(bed, min.gapwidth = 500, ignore.strand = FALSE)
#' collapseBed(bed, min.gapwidth = 500, ignore.strand = FALSE, return.idx.only = TRUE)
#'
#' @export
collapseBed <- function(bed,
                        min.gapwidth= 1L,
                        return.idx.only= FALSE,
                        ignore.strand= TRUE)
{
  # Import for incapsulation ----
  bed <- vlite::importBed(bed)

  # Reduce ----
  gr <- GenomicRanges::GRanges(bed)
  coll <- GenomicRanges::reduce(gr,
                                ignore.strand= ignore.strand,
                                min.gapwidth = min.gapwidth)

  if(return.idx.only)
  {
    # Return run-length id ----
    idx <- GenomicRanges::findOverlaps(gr,
                                       coll,
                                       ignore.strand= ignore.strand)
    idx <- subjectHits(idx)
    return(idx)
  } else {

    # Order and return collapsed ranges ----
    res <- vlite::importBed(coll)
    setorderv(res,
              c("seqnames", "start", "end"))
    return(res)
  }
}

#' Find Closest Genomic Features Between Two Sets of Regions
#'
#' @description
#' A wrapper around ?GenomicRanges::nearest and ?GenomicRanges::nearestKNeighbors that computes,
#' for each genomic range in a, the nearest genomic ranges in b and the stranded distance between them.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param k If specified, the k nearest neighbors are computed. If set to NA (default), only the closest
#' regions are returned, including ties.
#' @param ignore.strand If set to FALSE, only closest features that are on the same strand are reported.
#' If set to TRUE (default), closest features are reported on both strands. See details.
#'
#' @details
#' **Distance Calculation**:
#' - Features that have no nearest neighbor (e.g., chromsome missing in b) are not returned
#' - For overlapping features, a distance of 0 is returned.
#' - For non-overlapping features, the genomic distance between their closest boundaries are returned:
#'   * Negative distances indicate upstream features.
#'   * Positive distances indicate downstream features.
#'   * Unstranded features are treated as '+'.
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item idx.a: line index of the regions in a.
#'   \item idx.b: line index of the closest region(s) in b.
#'   \item dist: distance between closest regions (see details).
#' }
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr2R:100-200:+", "chr2R:100-200:-"))
#' b <- importBed(c("chr2R:300-400:+", "chr2R:300-400:-", "chr2R:3000-3100:+"))
#'
#' # Find single closest features
#' closestBed(a, b)[]
#'
#' Only consider features that are on the same strand
#' closestBed(a, b, ignore.strand = F)[]
#'
#' # Return all features at the second closest distance (including ties):
#' closestBed(a, b, k= 3)[]
#' closestBed(a, b, k= 3, ignore.strand= FALSE)[]
#'
#' @export
closestBed <- function(a,
                       b,
                       k= NA,
                       ignore.strand= TRUE)
{
  # Checks ----
  if(!is.na(k) && (k<1L | k %% 1 != 0))
    stop("If specified, k should be a positive integer")

  # Import a and b for encapsulation ----
  a <- vlite::importBed(a)
  b <- vlite::importBed(b)

  # Find nearest regions (n=1) ----
  if(is.na(k)) {
    closest <- suppressWarnings(
      GenomicRanges::nearest(
        GenomicRanges::GRanges(a),
        GenomicRanges::GRanges(b),
        select= "all",
        ignore.strand= ignore.strand
      )
    )
    closest <- as.data.table(closest)
    setnames(closest, c("idx.a", "idx.b"))
  } else {

    # Find k nearest neighbors ----
    a[, idx.a:= .I]
    b[, idx.b:= .I]
    # For each chromosome shared between a and b
    common.chr <- intersect(unique(a$seqnames), unique(b$seqnames))
    closest <- list()
    for(chr in common.chr) {
      .a <- a[seqnames==chr]
      .b <- b[seqnames==chr]
      # Find k nearest neighbors
      .cl <- GenomicRanges::nearestKNeighbors(
        GRanges(.a),
        GRanges(.b),
        k= k,
        ignore.strand= ignore.strand
      )
      .cl <- data.table(
        idx.a= rep(seq(.cl), lengths(.cl)),
        idx.b= unlist(.cl)
      )
      .cl[, idx.a:= .a$idx.a[idx.a]]
      .cl[, idx.b:= .b$idx.b[idx.b]]
      closest[[chr]] <- .cl
    }
    # Bind all chromosomes
    closest <- rbindlist(closest)
  }

  # Compute distances ----
  closest[, dist:= fcase(
    b$end[idx.b]<a$start[idx.a], b$end[idx.b]-a$start[idx.a],
    a$end[idx.a]<b$start[idx.b], b$start[idx.b]-a$end[idx.a],
    default = 0L
  )]

  # Upstream feature -> neg distance ----
  closest[, dist:= ifelse(a$strand[idx.a]=="-", -dist, dist)]

  # Return ----
  return(closest)
}

#' Calculate Feature Coverage Across Genomic Intervals
#'
#' @description
#' A wrapper around ?GenomicRanges::countOverlaps that computes, for each genomic range in a,
#' the number of overlapping regions in b.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param maxgap A single integer specifying the maximum gap allowed between 2 ranges for them to
#' be considered as overlapping. Default= -1.
#' @param minoverlap A single integer specifying the minimum overlap between 2 ranges for them to
#' be considered as overlapping. Default= 0.
#' @param ignore.strand If set to FALSE, only features that are on the same strand will be counted.
#' If set to TRUE (default), overlapping feature are counted irrespective of their strand.
#'
#' @return A numeric vector of length nrow(a) corresponding, for each region in a, to the number
#' of overlapping regions in b. 0 means no overlaps.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-300:+", "chr1:400-500:-", "chr2:100-200:+"))
#' b <- importBed(c("chr1:50-150:+", "chr1:200-250:+", "chr1:400-450:-", "chr1:400-450:+", "chr3:100-200:+"))
#'
#' # Count overlapping features
#' covBed(a, b)
#' covBed(a, b, ignore.strand = FALSE)
#'
#' @export
covBed <- function(a,
                   b,
                   maxgap= -1L,
                   minoverlap= 0L,
                   ignore.strand= TRUE)
{
  # Import for incapsulation ----
  a <- vlite::importBed(a)
  b <- vlite::importBed(b)

  # Compute coverage ----
  cov <- suppressWarnings(
    GenomicRanges::countOverlaps(
      query = GRanges(a),
      subject = GRanges(b),
      maxgap = maxgap,
      minoverlap = minoverlap,
      ignore.strand= ignore.strand
    )
  )

  # Return ----
  return(cov)
}

#' Find Overlapping Regions Between Two Sets of Genomic Intervals
#'
#' @description
#' A wrapper around ?GenomicRanges::findOverlaps() that identifies, for each genomic range in a,
#' the overlapping regions in b.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param maxgap A single integer specifying the maximum gap allowed between 2 ranges for them to
#' be considered as overlapping. Default= -1.
#' @param minoverlap A single integer specifying the minimum overlap between 2 ranges for them to
#' be considered as overlapping. Default= 0.
#' @param ignore.strand If set to FALSE, only reports overlaps between regions that are on the same strand.
#' If set to TRUE (default), reports overlaps on both strands.
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item idx.a: line index of the region in a.
#'   \item idx.b: line index of the overlapping region in b (NA if no overlap).
#'   \item overlap.start: overlap starting position (NA if no overlap).
#'   \item overlap.end: overlap ending position (NA if no overlap).
#'   \item overlap.width: overlap width (0 if no overlap).
#' }
#' Of note, regions in a with no overlap in b are not returned.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-200:+", "chr1:300-400:-", "chr2:100-200:+"))
#' b <- importBed(c("chr1:150-250:+", "chr1:350-450:+", "chr1:500-600:+"))
#'
#' # Find all overlaps
#' overlapBed(a, b)
#'
#' # Find strand-specific overlaps
#' overlapBed(a, b, ignore.strand = FALSE)
#'
#' @export
overlapBed <- function(a,
                       b,
                       maxgap= -1L,
                       minoverlap= 0L,
                       ignore.strand= TRUE)
{
  # Import bed for incapsulation ----
  a <- vlite::importBed(a)[, c("seqnames", "start", 'end', "strand")]
  b <- vlite::importBed(b)[, c("seqnames", "start", 'end', "strand")]

  # Overlaps ----
  ov <- suppressWarnings(
    GenomicRanges::findOverlaps(
      query = GenomicRanges::GRanges(a),
      subject = GenomicRanges::GRanges(b),
      maxgap = maxgap,
      minoverlap = minoverlap,
      ignore.strand = ignore.strand
    )
  )
  ov <- as.data.table(ov)
  setnames(ov, c("idx.a", "idx.b"))

  # Compute overlap.start, overlap.end, overlap.width
  coor <- cbind(a[ov$idx.a, .(start, end)],
                b[ov$idx.b, .(b.start= start, b.end= end)])
  coor[, overlap.start:= ifelse(start>b.start, start, b.start)]
  coor[, overlap.end:= ifelse(end<b.end, end, b.end)]
  coor[, overlap.width:= overlap.end-overlap.start+1]
  coor[is.na(overlap.width), overlap.width:= 0]

  # Final result ----
  res <- cbind(ov,
               coor[, .(overlap.start, overlap.end, overlap.width)])

  # Return ----
  return(res)
}

#' Subset overlapping or non-overlapping Regions
#'
#' @description
#' A wrapper around ?GenomicRanges::countOverlaps that subsets the genomic ranges in a that do
#' (or do not) overlap region(s) in b.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param maxgap A single integer specifying the maximum gap allowed between 2 ranges for them to
#' be considered as overlapping. Default= -1.
#' @param minoverlap A single integer specifying the minimum overlap between 2 ranges for them to
#' be considered as overlapping. Default= 0.
#' @param invert If set to TRUE, returns regions in a that have no overlap(s) in b.
#'   If set to FALSE (default), overlapping regions are returned.
#' @param ignore.strand If set to FALSE, only overlapping features that are on the same strand are reported.
#' If set to TRUE (default), overlaps will be computed irrespective of the strand.
#'
#' @return A data.table containing:
#' \itemize{
#'   \item When invert = FALSE: regions from a that overlap region(s) in b.
#'   \item When invert = TRUE: regions from a that DO NOT overlap any region(s) in b.
#' }
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-200:+", "chr1:300-400:-", "chr1:700-800:-", "chr2:500-600:+"))
#' b <- importBed(c("chr1:100-200:+", "chr1:350-450:+", "chr1:700-800:+"))
#'
#' # Find all overlapping regions
#' intersectBed(a, b)
#' intersectBed(a, b, minoverlap = 75L)
#' intersectBed(a, b, ignore.strand= FALSE)
#'
#' # Find non-overlapping regions
#' intersectBed(a, b, invert = TRUE)
#' intersectBed(a, b, invert = TRUE, ignore.strand= FALSE)
#'
#' @export
intersectBed <- function(a,
                         b,
                         maxgap= -1L,
                         minoverlap= 0L,
                         invert= FALSE,
                         ignore.strand= TRUE)
{
  # Import bed for incapsulation ----
  a <- vlite::importBed(a)
  b <- vlite::importBed(b)

  # Overlaps ----
  tmp <- rev(make.unique(c(names(a), "ov")))[1]
  assign(
    tmp,
    GenomicRanges::countOverlaps(
      query = GRanges(a),
      subject = GRanges(b),
      maxgap = maxgap,
      minoverlap = minoverlap,
      ignore.strand= ignore.strand
    )
  )

  # Non-intersecting indexes ----
  res <- if(invert) {
    a[get(tmp)==0]
  } else {
    a[get(tmp)>0]
  }

  # Return, preserving original order ----
  return(res)
}

#' Subtract Genomic Regions from Reference Intervals
#'
#' @description
#' A wrapper around ?GenomicRanges::setdiff that subtracts the regions in b to the genomic ranges in a.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param minoverlap A single integer specifying the minimum overlap with regions in for them to
#' be subtracted. Default= 1L.
#' @param ignore.strand If set to FALSE, only subtracts features that are on the same strand.
#' If set to TRUE (default), subtracts overlapping feature(s) regardless of their strand(s).
#'
#' @return A data.table containing the remaining portions of a, after subtracting the regions
#' defined in b.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:200-300:-", "chr1:100-500:-"))
#' b <- importBed(c("chr1:200-300:+", "chr1:400-450:-", "chr1:425-475:-"))
#'
#' # Basic example
#' subtractBed(a, b)
#'
#' # Only subtract regions with similar strand
#' subtractBed(a, b, ignore.strand= FALSE)
#'
#' @export
subtractBed <- function(a,
                        b,
                        minoverlap= 1L,
                        ignore.strand= TRUE)
{
  # Import for incapsulation ----
  a <- vlite::importBed(a)
  b <- vlite::importBed(b)

  # Subtract features ----
  sub <- suppressWarnings(
    GenomicRanges::subtract(
      GenomicRanges::GRanges(a),
      GenomicRanges::GRanges(b),
      minoverlap= minoverlap,
      ignore.strand= ignore.strand
    )
  )

  # Retrieve row indices from a ----
  tmp <- rev(make.unique(c(names(a), "idx.a")))[1]
  assign(
    tmp,
    rep(seq(nrow(a)), lengths(sub))
  )

  # Subtracted result ----
  res <- cbind(
    as.data.table(unlist(sub))[, !"width"],
    a[get(tmp), !c("seqnames", "start", "end", "strand")]
  )

  # Return ----
  return(res)
}

#' Clip Genomic Regions to Defined Limits
#'
#' @description
#' A wrapper around ?GenomicRanges::findOverlaps that clips the regions in a that extend beyond
#' the regions in b.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param ignore.strand If set to FALSE and strand column is provided, only clips boundaries
#' that are on the same strand. If set to TRUE (default), clips boundaries regardless of their strand.
#'
#' @return A data.table containing the remaining portions of a, after clipping them using the boundaries
#' defined in b.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-500:+"))
#' b <- importBed(c("chr1:200-300:+", "chr1:400-450:-", "chr1:425-475:-"))
#'
#' # Basic example
#' clipBed(a, b)[]
#'
#' # Strand-specific
#' clipBed(a, b, ignore.strand = FALSE)[]
#'
#' @export
clipBed <- function(a,
                    b,
                    ignore.strand= TRUE)
{
  # Import for incapsulation ----
  a <- vlite::importBed(a)
  b <- vlite::importBed(b)

  # Collapse regions in b ----
  coll <- collapseBed(b, ignore.strand= ignore.strand)

  # Compute overlaps with a ----
  tmp <- rev(make.unique(c(names(a), "ov")))[1]
  assign(
    tmp,
    vlite::overlapBed(
      a,
      coll,
      ignore.strand= ignore.strand
    )
  )

  # Resize overlaps ----
  res <- a[get(tmp)$idx.a]
  res[, start:= get(tmp)$overlap.start]
  res[, end:= get(tmp)$overlap.end]

  # Return ----
  return(res)
}

