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
# bam <- "db/bam/ORFtag/CpG1_input_rep2.2_mm10_collapsed.bam"
# bam.file <- "db/bam/ORFtag/CpG1_input_rep2.2_mm10.bam"

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

# If bam file is specified, add the coverage ----
if(length(args)==5) {
  # Check bam file exists
  bam.file <- args[5]
  if(!file.exists(bam.file))
    stop("Specified bam.file file could not be found to compute coverage.")
  # Save bed as temp file
  tmp <- tempfile(tmpdir = dirname(bed_file), fileext = ".tmp.bed")
  rtracklayer::export(.c[, .(seqnames, start, end= start, strand)],
                      con = tmp)
  # Compute coverage
  cmd <- paste("/software/2020/software/bedtools/2.27.1-foss-2018b/bin/bedtools coverage -s -a", tmp, "-b", bam.file)
  # Coverage is in column 7
  .c$score <- data.table::fread(cmd = cmd, header = FALSE)$V7
  # Check that it worked
  if(any(.c$score)==0)
    stop("Some scores were equal to 0, which happens when the job got OOM killed.
         Increase RAM and try again.")
  file.remove(tmp)
}

# Compute start ----
.c[strand=="+", start:= start-1]
.c[strand=="-", start:= start+width]
.c[, strand:= ifelse(strand=="+", "-", "+")] # Reverse strand (iPCR)
.c[, end:= start]

# Remove non unique insertions ----
.c <- unique(.c)
.c <- na.omit(.c)
setorderv(.c,
          c("seqnames", "start", "end"))

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




