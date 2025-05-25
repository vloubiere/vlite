#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% c(4, 5)) {
  stop("Please specify:\n
       [required] 1/ Collapsed bam file \n
       [required] 2/ Path to a gtf file containing non-first exons used for assignment \n
       [required] 3/ Output bed file (.bed) \n
       [required] 4/ Output assignment_table prefix \n
       [required] 5/ An optional bam file for which the coverage per insertion will be returned in the score column of the collapsed bed file \n")
}

# Load packges ----
suppressMessages(library(Rsamtools, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))
suppressMessages(library(data.table, warn.conflicts = FALSE))

# Tests ----
# bam <- "/groups/stark/vloubiere/projects/sebastian/db/bam/ORFtag/CpG1_input_rep2.2_mm10_collapsed.bam"
# bam.file <- "/groups/stark/vloubiere/projects/sebastian/db/bam/ORFtag/CpG1_input_rep2.2_mm10.bam"

# Variables ----
bam <- args[1]
exons <- args[2]
bed_file <- args[3]
assignment_table <- args[4]

# Import and parse bam file ----
params <- Rsamtools::ScanBamParam(what=c("rname", "strand", "pos", "qwidth"))
par.names <- c("seqnames", "strand", "start", "width")
.c <- Rsamtools::scanBam(bam, param = params)[[1]]
.c <- as.data.table(.c)
setnames(.c, par.names)

# Compute insertions' starting position, upstream of read1 ----
.c[strand=="+", start:= start] # Insertion on - strand! See below
.c[strand=="-", start:= start+width-1] # Insertion on + strand! See below
.c[, end:= start]

# Reverse strand (inverse PCR) ----
if("*" %in% .c$strand)
  stop("Error: reads with strand '*' were detected, which can't be assigned.")
.c[, strand:= ifelse(strand=="+", "-", "+")]

# Remove non unique insertions (there shouldn't be any) ----
.c <- na.omit(unique(.c))
setorderv(.c,
          c("seqnames", "start", "end"))

# If bam file is specified, add the coverage ----
if(length(args)==5) {
  bam.file <- args[5]
  # Import full bam file
  bam.file <- Rsamtools::scanBam(bam.file, param = params)[[1]]
  bam.file <- as.data.table(bam.file)
  setnames(bam.file, par.names)
  # Compute start and reverse strand (similarly to input bam, see above)
  bam.file[strand=="+", start:= start]
  bam.file[strand=="-", start:= start+width-1]
  bam.file[, strand:= ifelse(strand=="+", "-", "+")]
  # Compute duplicate counts
  .c$score <- bam.file[.c, .N, .EACHI, on= c("seqnames", "start", "strand")]$N
}

# Save bed file ----
.c <- GenomicRanges::GRanges(.c)
rtracklayer::export(.c, bed_file)

# Import non-first exons gtf file ----
exons <- rtracklayer::import(exons)
cols <- c("gene_id", "gene_name", "mgi_id", "exon_number", "exon_id")
cols <- cols[cols %in% names(mcols(exons))] # Human gtf does not contain mgi_id
mcols(exons) <- mcols(exons[, cols])

# For each strand, assign insertions to closest downstream, non-first exon ----
for(.strand in c("same_strand", "rev_strand"))
{
  curr <- data.table::copy(.c)
  if("score" %in% names(curr))
    setnames(curr, "score", "ins_cov")
  if(.strand=="rev_strand")
    strand(curr) <- sapply(strand(curr), switch,  "-"="+", "+"="-")
  idx <- IRanges::precede(curr, exons, ignore.strand= FALSE)
  mcols(curr) <- NULL
  mcols(curr) <- mcols(exons)[idx,]
  mcols(curr)$dist <- NA
  curr$dist[!is.na(idx)] <- IRanges::distance(curr[!is.na(idx)], exons[idx[!is.na(idx)]])
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




