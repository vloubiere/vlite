#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% 4:5) {  # Fix the condition to match required arguments
  stop("Please specify:\n",
       "[required] 1/ Input bam file\n",
       "[required] 2/ Layout ('SINGLE' or 'PAIRED')\n",
       "[required] 3/ Should the fragments be extended? (TRUE/FALSE)\n",
       "[required] 4/ bw output file path\n",
       "[optional] 5/ Integer. Number of nt by which the reads should be extended before computing coverage.")
}

# Load packages ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(Rsamtools, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))

# Parse arguments ----
bam <- args[1]
layout <- args[2]
extend.fragment <- args[3]
bw <- args[4]
extsize <- ifelse(length(args)==5, args[5], 0)

# Input validation
stopifnot(file.exists(bam))
stopifnot(layout %in% c("SINGLE", "PAIRED"))
extend.fragment <- as.logical(extend.fragment)  # Convert string to logical
stopifnot(!is.na(extend.fragment))
extsize <- as.numeric(extsize)
stopifnot(extsize >= 0)

# Import bam file ----
param <- Rsamtools::ScanBamParam(what= c("qname", "flag", "rname", "pos", "qwidth", "strand"))
reads <- Rsamtools::scanBam(bam,
                            param = param)[[1]]
reads <- as.data.table(reads)
setnames(reads,
         c("qname", "flag", "rname", "pos", "qwidth", "strand"),
         c("readID", "flag", "seqnames", "start", "width", "strand"))
rm(param)
gc()

# Identify unambiguously mapped paired-end or single-end reads ----
if(layout=="PAIRED") {
  reads[, check:= (bitwAnd(flag, 0x4) == 0) & (bitwAnd(flag, 0x2) != 0)]
} else if(layout=="SINGLE") {
  reads[, check:=
          (bitwAnd(flag, 0x4) == 0) &
          (bitwAnd(flag, 0x100) == 0) &
          (bitwAnd(flag, 0x800) == 0)]
}
reads <- reads[(check)]
reads$check <- NULL

# Extend reads ----
reads[, end:= start+width-1]
if(layout=="PAIRED" && extend.fragment) {
  # Make sure the two mates are on the same chromosome
  reads[, sameChr:= .N == 2, .(seqnames, readID)]
  # Extend fragment
  reads <- reads[(sameChr), .(start= min(start), end= max(end)), .(seqnames, readID)]
} else if(extsize>0) {
  # Add extsize depending on strand
  reads[strand!="-", end:= end+extsize]
  reads[strand=="-", start:= start-extsize]
}

# Get chromosome sizes from BAM header for proper bigwig format ----
chr_sizes <- Rsamtools::scanBamHeader(bam)[[1]]$targets

# Clip to chromosome limits ----
reads[as.data.table(chr_sizes, keep.rownames = "seqnames"), max:= i.chr_sizes, on= "seqnames"]
reads[start<1, start:= 1]
reads[end>max, end:= max]

# Create GRange ----
gr <- GenomicRanges::GRanges(
  seqnames = reads$seqnames,
  ranges = IRanges::IRanges(start = reads$start, end = reads$end),
  seqlengths = chr_sizes
)

# Compute coverage and export ----
cov <- GenomicRanges::coverage(gr)
rtracklayer::export(cov, bw)
