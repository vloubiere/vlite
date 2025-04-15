#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript script.R input_file.bdg output_file.bw genome")
}

# Load packages ----
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomeInfoDb, warn.conflicts = FALSE))
suppressMessages(library(BSgenome, warn.conflicts = FALSE))

# Parse arguments ----
bdg.file <- args[1]
bw.file <- args[2]
genome <- args[3]

# Import bedgraph
gr <- rtracklayer::import(bdg.file, format = "bedGraph")
chrLengths <- GenomeInfoDb::seqlengths(BSgenome::getBSgenome(genome))
seqlengths(gr) <- chrLengths[seqlevels(gr)]

# Export to .bw file
rtracklayer::export(gr, bw.file, format = "bigWig")
