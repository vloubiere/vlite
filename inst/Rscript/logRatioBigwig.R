#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 8) {  # Fix the condition to match required arguments
  stop("Please specify:\n",
       "[required] 1/ Bed file containing experiment reads. \n",
       "[required] 2/ Bed file containing input/control reads. \n",
       "[required] 3/ A BS genome name. \n",
       "[required] 4/ Optional path to a bed file. If provided, only reads overlapping these regions on the same strand are used. \n",
       "[required] 5/ Should a gaussian smoothing be applied? \n",
       "[required] 6/ Pseudocount. If set to NULL empty, it will be set to the 0.01 percentile of non-0 value \n",
       "[required] 7/ Output .bw file path \n",
       "[required] 8/ Bins width \n")
}

# Load packges ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomeInfoDb, warn.conflicts = FALSE))
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Parse arguments
track <- args[1]
input.track <- args[2]
genome <- args[3]
regions <- if(args[4] == "NULL") NULL else args[4] # ifelse does not work here!
gaussian <- as.logical(args[5])
pseudocount <- if(args[6] == "NULL") NULL else as.numeric(args[5]) # ifelse does not work here!
output.file <- args[7]
bins.width <- as.integer(args[8])

# Checks
stopifnot(grepl(".bed$", track))
stopifnot(file.exists(track))
stopifnot(grepl(".bed$", input.track))
stopifnot(file.exists(input.track))
if(!is.null(regions)) {
  stopifnot(grepl(".bed$", regions))
  stopifnot(file.exists(regions))
}

# Import bed files ----
signal <- rtracklayer::import(track)
input <- rtracklayer::import(input.track)

# Restrict to regions ----
if(!is.null(regions)) {
  regions <- rtracklayer::import(regions)
  signal <- signal[countOverlaps(signal, regions, ignore.strand= FALSE)>0]
  seqlevels(signal) <- as.character(unique(seqnames(signal)))
  input <- input[countOverlaps(input, regions, ignore.strand= FALSE)>0]
  seqlevels(input) <- as.character(unique(seqnames(input)))
}

# Bins covered regions ----
cov.regions <- GenomicRanges::reduce(c(signal, input), ignore.strand= TRUE)
bins <- GenomicRanges::slidingWindows(cov.regions, width = bins.width, step = bins.width)
bins <- unlist(bins)

# Count Overlaps ----
bins$signal <- GenomicRanges::countOverlaps(bins, signal, ignore.strand= TRUE) / length(signal) * 1e6
bins$input <- GenomicRanges::countOverlaps(bins, input, ignore.strand= TRUE) / length(input) * 1e6

# Gaussian smoothing ----
if(gaussian) {
  bins$signal <- vlite::gaussianBlur(bins$signal, size = 5, sigma = 1)
  bins$input <- vlite::gaussianBlur(bins$input, size = 5, sigma = 1)
}

# Add pseudocount if 0s in the data ----
if(any(c(bins$signal, bins$input)==0)) {
  if(is.null(pseudocount)) {
    pseudocount <- c(bins$signal, bins$input)
    pseudocount <- quantile(pseudocount[pseudocount>0], 0.01)
  }
  bins$signal <- bins$signal+pseudocount
  bins$input <- bins$input+pseudocount
}

# Calculate log2 ratio
bins$score <- log2(bins$signal / bins$input)

# Add seqLengths ----
chrLengths <- GenomeInfoDb::seqlengths(BSgenome::getBSgenome(genome))
seqlengths(bins) <- chrLengths[seqlevels(bins)]

# Export as BigWig
rtracklayer::export(bins, output.file, format = "BigWig")
