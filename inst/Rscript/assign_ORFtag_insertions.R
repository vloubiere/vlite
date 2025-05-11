#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) == 4) {
  stop("Please specify:\n
       [required] 1/ Collapsed bam file\n
       [required] 2/ Path to a gtf file containing non-first exons used for assignment\n
       [required] 3/ Output bed file (.bed)\n
       [required] 4/ Output assignment_table prefix\n")
}

# Load packges ----
suppressMessages(library(Rsamtools, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))
suppressMessages(library(data.table, warn.conflicts = FALSE))

# Variables ----
bam <- args[1]
exons <- args[2]
bed_file <- args[3]
assignment_table <- args[4]

# Import and parse bam file ----
.c <- Rsamtools::scanBam(bam,
                         param = Rsamtools::ScanBamParam(what=c("rname", "strand", "pos", "qwidth")))[[1]]
.c <- as.data.table(.c)
setnames(.c, c("seqnames", "strand", "start", "width"))
# Compute start
.c[strand=="+", start:= start-1]
.c[strand=="-", start:= start+width]
.c[, strand:= ifelse(strand=="+", "-", "+")]# Reverse strand (iPCR)
.c[, end:= start]
# Remove non unique insertions
.c <- unique(.c)
.c <- na.omit(.c)
setorderv(.c, c("seqnames", "start", "end"))
# Save bed file
.c <- GenomicRanges::GRanges(.c)
rtracklayer::export(.c, bed_file)
# Import non-first exons gtf file
exons <- rtracklayer::import(exons)
cols <- c("gene_id", "gene_name", "mgi_id", "exon_number", "exon_id")
cols <- cols[cols %in% names(mcols(exons))] # Human gtf does not contain mgi_id
mcols(exons) <- mcols(exons[, cols])
# For each strand, assign insertions to closest downstream, non-first exon ----
for(.strand in c("same_strand", "rev_strand"))
{
  curr <- .c
  if(.strand=="rev_strand")
    strand(curr) <- sapply(strand(curr), switch,  "-"="+", "+"="-", "*"="*")
  idx <- precede(curr, exons)
  mcols(curr) <- NULL
  mcols(curr) <- mcols(exons)[idx,]
  mcols(curr)$dist <- NA
  curr$dist[!is.na(idx)] <- distance(curr[!is.na(idx)], exons[idx[!is.na(idx)]])
  # SAVE
  curr <- as.data.table(curr)
  curr$width <- NULL
  fwrite(curr,
         paste0(assignment_table, "_", .strand, ".txt"),
         col.names = T,
         quote= F,
         sep= "\t",
         row.names = F,
         na= NA)
}




