#' Get BS genome sizes
#'
#' @param genome BSgenome object to use: "dm6", "mm10"...
#' @examples
#' getBSgenomeSize("mm10")
#'
#' @return data.table containing chromosomes start and end coordinates
#' @export
getBSgenomeSize <- function(genome)
{
  genome <- BSgenome::getBSgenome(genome, load.only = TRUE)
  gSize <- GenomeInfoDb::seqinfo(genome)
  gSize <- data.table(seqnames= GenomicRanges::seqnames(gSize),
                      start= 1L,
                      end= GenomeInfoDb::seqlengths(gSize))
  gSize[, strand:= "*"]
  return(gSize)
}

#' Randomly Sample Control Genomic Regions
#'
#' @description
#' Generates random genomic regions from a BSgenome, matching specified widths or the distribution of an input BED file.
#'
#' @param genome BSgenome name: "dm6", "mm10"...
#' @param widths Integer vector specifying the width of the regions to sample. The length of the vector determines the number of regions.
#' @param bed Optional genomic ranges. If provided, the widths arguments will be ignored, and control regions will match be sampled to match
#' the number of regions in bed, as well as their widths and seqnames distribution.
#' @param restrict.seqnames Optional character vector. Restrict sampling to these seqnames.
#' @param no.overlap If set to TRUE (default), ensures that there is no overlap between sampled and bed regions (if provided).
#' Stops after 10 iterations.
#'
#' @return A gr data.table of sampled control regions.
#'
#' @examples
#' # Sample using widths
#' test <- randomRegionsBSgenome("dm6", widths= rep(5000, 1e3))
#'
#' Sampling using bed
#' randomRegionsBSgenome(genome = "dm6", bed = test)
#'
#' @export
randomRegionsBSgenome <- function(genome,
                                  widths= rep(100, 1000),
                                  bed,
                                  restrict.seqnames= NULL,
                                  no.overlap= TRUE)
{
  # Checks ----
  if(!is.null(restrict.seqnames) & !missing(bed))
    stop("When bed is provided, restrict.seqnames should be set to NULL.")
  if(!missing(bed)) {
    bed <- importBed(bed)
    widths <- bed[, end-start+1]
  }

  # Retrieve chromosome sizes ----
  gSize <- getBSgenomeSize(genome= genome)

  # Remove chromosomes that are too short ----
  too.short <- gSize$end < max(widths)
  if(any(too.short)) {
    message(paste(sum(too.short), "/", nrow(gSize), "chromosomes shorter than max.width -> removed!"))
    gSize <- gSize[!(too.short)]
  }

  # Restrict seqnames ----
  if(!is.null(restrict.seqnames))
    gSize <- gSize[seqnames %in% restrict.seqnames]

  # Compute chromosome probabilities ----
  prob <- if(missing(bed)) gSize$end else gSize[, sum(bed$seqnames==seqnames), seqnames]$V1

  # Random sampling ----
  rdm <- rdmSamplingBS(gSize= gSize,
                       prob= prob,
                       widths= widths)

  # Avoid overlaps ----
  if(no.overlap && !missing(bed)) {
    remove <- covBed(rdm, bed)>0
    nTries <- 1
    while(any(remove) & nTries<10) {
      add <- rdmSamplingBS(gSize= gSize,
                           prob= prob,
                           widths= rdm$width[remove])
      rdm <- rbind(rdm[!(remove)], add)
      remove <- covBed(rdm, bed)>0
      remaining <- sum(remove)
      nTries <- nTries+1
    }
    if(remaining>0)
      warning(paste(remaining, "overlaps could not be fixed after 10 attemps and will me missing"))
    rdm <- intersectBed(rdm, bed, invert = T)
  }

  # Return result ----
  rdm <- rdm[, .(seqnames, start, end, strand)]
  return(rdm)
}

#' Get genomic sequence
#'
#' Returns the sequences of a bed.
#'
#' @param bed Input regions in any format compatible with ?importBed.
#' @param genome BSgenome name: "dm6", "mm10"...
#'
#' @examples
#' # Sample random regions
#' test <- randomRegionsBSgenome("dm6", widths= rep(100, 3))
#' getBSsequence(test, "dm6")
#'
#' @return Character vector of sequences
#' @export
getBSsequence <- function(bed,
                          genome)
{
  # Import ----
  bed <- importBed(bed)

  # Make sequence names ----
  if(!"name" %in% names(bed)) {
    bed[, name:= paste0(seqnames, ":", start, "-", end, ":", strand)]
  }

  # Extract sequences ----
  sequences <- BSgenome::getSeq(BSgenome::getBSgenome(genome, load.only = TRUE),
                                names= bed$seqnames,
                                start= bed$start,
                                end= bed$end,
                                strand= bed$strand,
                                as.character= T)

  # Add names and return ----
  names(sequences) <- bed$name
  return(sequences)
}
