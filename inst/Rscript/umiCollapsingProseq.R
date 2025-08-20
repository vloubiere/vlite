#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("Please specify:\n
       [required] 1/ Aligned bam file containing UMIs\n
       [required] 2/ Output count file (.txt)\n
       [required] 3/ Output stats file (.txt)\n
       [required] 4/ Should the strand of the read be flipped? Default= TRUE for PRO-Seq\n")
}

suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(Rsamtools, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))
suppressMessages(library(stringdist, warn.conflicts = FALSE))

# Variables ----
bam <- args[1]
counts.output <- args[2]
stats.output <- args[3]
flip.strand <- as.logical(args[4])

# Import data ----
param <- Rsamtools::ScanBamParam(what= c("qname", "rname", "pos", "strand", "qwidth", "mapq"))
dat <- Rsamtools::scanBam(bam,
                          param = param)[[1]]
dat <- as.data.table(dat)
setnames(dat,
         c("qname", "rname", "pos", "strand", "qwidth", "mapq"),
         c("read", "seqnames", "start", "strand", "width", "mapq"))

# Resize to single nucleotide ----
dat[, strand:= as.character(strand)]
dat[strand=="-", start:= start+width-1]

# Flip the strand of the reads ----
if(flip.strand) # For PRO-Seq, should be set to TRUE. FALSE for STAP-Seq
  dat[!is.na(strand), strand:= ifelse(strand=="+", "-", "+")]

# Save coordinates ----
dat[, coor:= paste0(seqnames, ":", start, ":", strand)]

# Extract UMIs ----
dat[!is.na(seqnames), UMI:= sub(".*_(.*)$", "\\1", read)]
dat[is.na(seqnames), UMI:= "GGGGGGGGGG"] # Unaligned reads (not collapsed, keep for total reads)
umi_length <- nchar(dat$UMI)
if(any(umi_length<10)) {
  stop("Some UMIs are shorter than 10nt")
} else if(any(umi_length>10)) {
  warning("Some UMIs were longer than 10nt and will be trimmed.")
  dat[umi_length>10, UMI:= substr(UMI, 1, 10)]
}

# Compute total counts ----
dat <- dat[, .(umi_N= .N), .(coor, UMI)]
dat[, total_counts:= sum(umi_N), coor]
setorderv(dat, "umi_N", order = -1)

# Check whether UMI might be collapsed ----
dat[, collapsed:= T, coor]
dat[, idx:= .I]
for(i in 1:10)
{
  dat[, check:= idx[1], .(coor, gsub(paste0("^(.{", i-1, "})."), "\\1", UMI))]
  potentialDup <- unique(dat[(check<idx), c(check, idx)])
  dat[potentialDup, collapsed:= FALSE]
  print(i)
}
dat$idx <- NULL
paste0(sum(dat$collapsed), " / ", nrow(dat), " pre-collapsed")

# UMI collapsing (>1 diff) ----
while(any(!dat$collapsed))
{
  dat[!(collapsed), c("collapsed", "UMI"):= {
    coll <- stringdist(UMI[1],
                       UMI,
                       method="hamming",
                       nthread= getDTthreads()-1)<=1
    UMI[coll] <- UMI[1]
    .(coll, UMI)
  }, coor]
}

# Final collapsing ----
dat <- unique(dat[, .(coor, total_counts, UMI)])
dat <- dat[, .(umi_counts= .N), .(coor, total_counts)]
dat[coor=="NA:NA:NA", umi_counts:= NA]

# Save output ----
fwrite(dat,
       counts.output,
       sep= "\t",
       quote= F,
       na= NA)

# Compute statistics and spikein sizeFactor ----
stats <- data.table(total= sum(dat$total_counts),
                    mapped= sum(dat[coor!="NA:NA:NA", total_counts]), # Remove unmapped
                    umi_counts= sum(dat[coor!="NA:NA:NA", umi_counts]))
fwrite(stats,
       stats.output,
       sep= "\t",
       na= NA)
