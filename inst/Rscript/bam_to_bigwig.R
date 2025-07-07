#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% 5:6) {  # Fix the condition to match required arguments
  stop("Please specify:\n",
       "[required] 1/ Input bam file \n",
       "[required] 2/ Layout ('SINGLE' or 'PAIRED') \n",
       "[required] 3/ Should the output bigwig be CPM normalized? \n",
       "[required] 4/ Should the fragments be extended? (TRUE/FALSE) \n",
       "[required] 5/ bw output file path \n",
       "[optional] 6/ Integer. Number of nt by which the reads should be extended before computing coverage.")
}

# Load packages ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(Rsamtools, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))

# Parse arguments ----
bam <- args[1]
layout <- args[2]
libsize.norm <- as.logical(args[3])
extend.fragment <- as.logical(args[4])
bw <- args[5]
extsize <- ifelse(length(args)==6, args[6], 0)

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

# Filter mapped paired-end or single-end reads ----
if (layout == "PAIRED") {
  reads[, check := (bitwAnd(flag, 0x4) == 0) &    # mapped
                   (bitwAnd(flag, 0x2) == 0x2) &  # properly paired
                   (bitwAnd(flag, 0x100) == 0) &  # not secondary
                   (bitwAnd(flag, 0x800) == 0)]   # not supplementary
} else if (layout == "SINGLE") {
  reads[, check := (bitwAnd(flag, 0x4) == 0) &   # mapped
                   (bitwAnd(flag, 0x100) == 0) & # not secondary
                   (bitwAnd(flag, 0x800) == 0)]  # not supplementary
}
reads <- reads[(check)]
reads$check <- NULL
setkeyv(reads, "readID")

# Extend reads ----
reads[, end:= start+width-1]
if(layout=="PAIRED" && extend.fragment) {
  # Make sure the two mates are on the same chromosome
  reads[, Nchr:= .N, .(seqnames, readID)]
  # Extend fragment
  reads <- reads[(Nchr==2), .(start= min(start), end= max(end)), .(seqnames, readID)]
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

# Create GRanges ----
gr <- GenomicRanges::GRanges(
  seqnames = reads$seqnames,
  ranges = IRanges::IRanges(start = reads$start, end = reads$end),
  seqlengths = chr_sizes
)

# Compute coverage and normalize ----
cov <- GenomicRanges::coverage(gr)
if(libsize.norm)
  cov <- cov/nrow(reads)*1e6

# Export ----
rtracklayer::export(cov, bw)
