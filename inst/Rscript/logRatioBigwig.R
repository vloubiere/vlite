#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly = TRUE)
if(!length(args) == 5) {  # Fix the condition to match required arguments
  stop("Please specify:\n",
       "[required] 1/ Experiment bw file. \n",
       "[required] 2/ Input bw file. \n",
       "[required] 3/ An optional bed file to restict regions for which log2Ratio will be computed \n",
       "[required] 4/ Pseudocount. If left empty, it will be set to 1% percentile of non-0 values \n",
       "[required] 5/ Output .bw file path \n")
}

# Load packges ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomeInfoDb, warn.conflicts = FALSE))
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Parse arguments
track <- args[1]
input <- args[2]
bed <- if(args[3] == "NULL") NULL else args[3] # ifelse does not work here!
pseudocount <- if(args[4] == "NULL") NULL else as.numeric(args[4]) # ifelse does not work here!
output.file <- args[5]

# Restric to loci of interest
if(!is.null(bed))
  regions <- rtracklayer::import(bed)

# Import experiment signal)
signal <- if(!is.null(bed))
  rtracklayer::import(track, which = regions, as = "RleList") else
    rtracklayer::import(track, as = "RleList")

# Import input signal
input <- if(!is.null(bed))
  rtracklayer::import(input, which = regions, as = "RleList") else
    rtracklayer::import(track, as = "RleList")

# Restrict to common chromosomes
common_chroms <- intersect(names(signal), names(input))
signal <- signal[common_chroms]
input <- input[common_chroms]

# Log2 ratio
if(is.null(pseudocount))
  pseudocount <- min(sapply(c(input, signal), function(x) quantile(x[x>0], 0.01)))
log2_ratio <- log2((signal+pseudocount) / (input+pseudocount))

# Save
rtracklayer::export(log2_ratio, output.file, format = "BigWig")
