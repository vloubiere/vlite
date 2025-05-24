#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) {  # Fix the condition to match required arguments
  stop("Please specify:\n",
       "[required] 1/ Input bam file, where the UMI sequence should be appended to the read ID\n",
       "[required] 2/ Layout ('SINGLE' or 'PAIRED')\n",
       "[required] 3/ umi counts output file (.txt)\n",
       "[required] 4/ umi-collapsed bed output (.bed)\n")
}

# Load packages ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(Rsamtools, warn.conflicts = FALSE))
suppressMessages(library(stringdist, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))

# Parse arguments ----
# bam <- "/groups/stark/vloubiere/projects/Facilitators/db/bam/screen_250bp_K562_rep1_mm9.bam"
bam <- args[1]
layout <- args[2]
counts.output.file <- args[3]
bed.output.file <- args[4]

# Input validation
stopifnot(file.exists(bam))
stopifnot(layout %in% c("SINGLE", "PAIRED"))
stopifnot(grepl(".txt$", counts.output.file))
stopifnot(grepl(".bed$", bed.output.file))

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
  reads[, check := (bitwAnd(flag, 0x4) == 0) &  # mapped
          (bitwAnd(flag, 0x2) == 0x2)] # properly paired
} else if (layout == "SINGLE") {
  reads[, check := (bitwAnd(flag, 0x4) == 0)] # mapped
}
reads <- reads[(check)]
reads$check <- NULL
setkeyv(reads, "readID")

# Extend reads ----
reads[, end:= start+width-1]
if(layout=="PAIRED") {
  # Make sure the two mates are on the same chromosome
  reads[, Nchr:= .N, .(seqnames, readID)]
  # Retrieve first mate per pair
  reads[, isFirst:= bitwAnd(flag, 0x40) == 0x40, flag]
  # Extend fragment
  reads <- reads[(Nchr==2), .(start= min(start), end= max(end), strand= strand[isFirst]), .(seqnames, readID)]
}

# Simplify coordinates ----
reads[, coor:= paste0(seqnames, ":", start, "-", end, ":", strand)]

# Check that UMIs have been appended to the readID ----
if(any(!grepl("_[ATCGN]+$", reads$readID)))
  stop("Some read IDs have no UMI sequences attached.")

# Extract UMIs ----
reads[, UMI:= gsub(".*_([A-Z]{10}).*", "\\1", readID)]
reads <- reads[, .(umi_N= .N), .(coor, UMI)]
reads[, total_counts:= sum(umi_N), .(coor)]
setorderv(reads, "umi_N", order = -1)

# Quick check to identify UMIs that might require collapsing
reads[, collapsed:= T, coor]
reads[, idx:= .I]
for(i in 1:10)
{
  reads[, check:= idx[1], .(coor, gsub(paste0("^(.{", i-1, "})."), "\\1", UMI))]
  potentialDup <- unique(reads[(check<idx), c(check, idx)])
  reads[potentialDup, collapsed:= FALSE]
  print(i)
}
reads$idx <- NULL
paste0(sum(reads$collapsed), " / ", nrow(reads), " pre-collapsed")

# UMI collapsing (>1 diff) ----
while(any(!reads$collapsed))
{
  reads[!(collapsed), c("collapsed", "UMI"):= {
    coll <- stringdist::stringdist(
      UMI[1],
      UMI,
      method="hamming",
      nthread= getDTthreads()-1
    )<=1
    UMI[coll] <- UMI[1]
    .(coll, UMI)
  }, coor]
}

# Final collapsing ----
final <- unique(reads[, .(coor, total_counts, UMI)])
final <- final[, .(umi_counts= .N), .(coor, total_counts)]

# SAVE ----
fwrite(x = final,
       file = counts.output.file,
       sep= "\t",
       quote= FALSE,
       row.names = FALSE,
       col.names = TRUE)

# Bed ----
bed <- final[rep(seq(nrow(final)), final$umi_counts)]
bed[, c("seqnames", "coor", "strand"):= tstrsplit(coor, ":", type.convert = TRUE)]
bed[, c("start", "end"):= tstrsplit(coor, "-", type.convert = TRUE)]
bed <- bed[, .(seqnames, start, end, strand)]
setorderv(bed,
          c("seqnames", "start", "end", "strand"))

# SAVE ----
rtracklayer::export(object = bed,
                    con = bed.output.file)
