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

#' Randomly Sample Control Genomic Regions
#'
#' @description
#' Randomly sample genomic regions from a BSgenome, matching specified widths and seqnames.
#'
#' @param genome BSgenome name: "dm6", "mm10"...
#' @param widths Integer vector specifying the width of the regions to sample. The length of the vector
#' determines the number of sampled regions.
#' @param restrict.seqnames If specified, restricts random sampling to the corresponding seqnames.
#' @param no.overlaps Genomic ranges for which overlaps should be avoided dhring sampling, in any format compatible
#' with ?importBed. Default= NULL.
#' @param ignore.strand If genomic ranges are provided in no.overlaps, only overlaps on the same strand will be
#' avoided. If set to TRUE (default), overlaps on both strands are avoided.
#'
#' @return A gr data.table of sampled control regions.
#'
#' @examples
#' # Sample using widths
#' test <- randomRegionsBSgenome(genome= "dm6",
#' widths= sample(c(1000, 2000, 5000), 1e3, replace= T),
#' restrict.seqnames= c("chr2L", "chr2R", "chr3L", "chr3R"))
#'
#' Sampling using bed
#' rdm <- randomRegionsBSgenome(genome = "dm6", no.overlaps = test)
#'
#' @export
randomRegionsBSgenome <- function(genome,
                                  widths,
                                  restrict.seqnames= NULL,
                                  no.overlaps= NULL,
                                  ignore.strand= TRUE)
{
  # Checks ----
  if(!is.null(restrict.seqnames) && length(restrict.seqnames)==0)
    stop("restrict.seqnames should contain at least one seqname.")

  # Retrieve chromosome sizes ----
  bed <- getBSgenomeSize(genome= genome)

  # Restrict seqnames ----
  if(!is.null(restrict.seqnames)) {
    if(!all(restrict.seqnames %in% bed$seqnames))
      stop("All provided restrict.seqnames should exist in the provided genome.")
    bed <- bed[seqnames %in% restrict.seqnames]
  }

  # Random sampling ----
  rdm <- randomRegionsBed(bed= bed,
                          widths= widths,
                          no.overlaps= no.overlaps,
                          ignore.strand= ignore.strand)

  # Return ----
  return(rdm)
}
