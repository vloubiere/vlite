#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Please specify:\n",
       "[required] 1/ Input bedgraph (.bdg) file \n",
       "[optional] 2/ Output bigwig (.bw) file \n",
       "[required] 3/ A BS genome name ('dm6', 'mm11', ...) \n",
       "[required] 4/ A scaling factor by which the bedgraph signal will be divided \n")
}

# Load packages ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomeInfoDb, warn.conflicts = FALSE))
suppressMessages(library(BSgenome, warn.conflicts = FALSE))

# Parse arguments ----
bdg.file <- args[1]
bw.file <- args[2]
genome <- args[3]
scaling.factor <- as.numeric(args[4])

# Import bedgraph ----
gr <- rtracklayer::import(bdg.file, format = "bedGraph")

# Apply scaling factor ----
if(scaling.factor!=1)
  gr$score <- gr$score/scaling.factor

# Add seqLengths ----
chrLengths <- GenomeInfoDb::seqlengths(BSgenome::getBSgenome(genome))
seqlengths(gr) <- chrLengths[seqlevels(gr)]

# Export as BigWig
rtracklayer::export(gr, bw.file, format = "BigWig")
