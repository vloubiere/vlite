#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) == 3) {  # Fix the condition to match required arguments
  stop("Please specify:\n",
       "[required] 1/ Input bed file\n",
       "[optional] 2/ Genome\n",
       "[required] 3/ bw output file path\n")
}

# Load packages ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))

# Parse arguments ----
bed <- args[1]
genome <- args[2]
bw <- args[3]

# Input validation
stopifnot(file.exists(bed))
stopifnot(grepl(".bw$", bw))

# Import bed file ----
gr <- rtracklayer::import(bed)

# Add genome ----
chrLengths <- GenomeInfoDb::seqlengths(BSgenome::getBSgenome(genome))
seqlengths(gr) <- chrLengths[seqlevels(gr)]

# Compute coverage ----
norm_cov <- GenomicRanges::coverage(gr) / length(gr) * 1e6

# Compute coverage and export ----
rtracklayer::export(norm_cov, bw)
