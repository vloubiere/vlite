#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2) {
  stop("Please specify:\n
       [required] 1/ Aligned BCs bam\n
       [required] 3/ Output counts file (.txt)\n")
}

# Load packges ----
suppressMessages(library(Rsamtools, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))
suppressMessages(library(data.table, warn.conflicts = FALSE))

# Variables ----
bam <- args[1]
output_counts <- args[2]

# Import bam ----
.c <- Rsamtools::scanBam(bam,
                         param = Rsamtools::ScanBamParam(what= c("qname", "rname", "mapq", "strand", "seq")))
.c <- .c[[1]]
.c$seq <- as.character(.c$seq)
.c <- as.data.table(.c)
.c[, rname:= as.character(rname)]
.c[, strand:= as.character(strand)]

# Count ----
.c <- .c[, .(count= .N), rname]

# Order based on alignment index ----
hdr <- Rsamtools::scanBamHeader(bam)
dic <- data.table(seqnames = names(hdr[[1]]$targets))
dic[.c, count:= i.count, on= "seqnames==rname"]
dic[is.na(count), count:= 0]

# Save ----
fwrite(dic,
       output_counts,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
