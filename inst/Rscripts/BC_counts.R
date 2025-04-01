#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) == 3) {
  stop("Please specify:\n
       [required] 1/ Aligned BCs bam\n
       [required] 2/ Dictionary file (.rds)\n
       [required] 3/ Output counts file (.txt)\n")
}

# Load packges ----
suppressMessages(library(Rsamtools, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))
suppressMessages(library(data.table, warn.conflicts = FALSE))

# Variables ----
bam <- args[1]
BC <- args[2]
output_counts <- args[3]

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

# Import dictionary ----
dic <- readRDS(BC)
setnames(dic, "bcID", "ID")
dic <- unique(dic[, .(ID, BC)])
dic[.c, count:= i.count, on= "ID==rname"]
dic[is.na(count), count:= 0]

# Save ----
fwrite(dic,
       output_counts,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
