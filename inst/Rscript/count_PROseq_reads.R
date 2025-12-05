#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there 3-4 args are returned. If not, return an error
if (!length(args) %in% c(3, 4)) {
  stop("Please specify:\n
       [required] 1/ Reference genome umi count file \n
       [required] 2/ A .rds annotation file \n
       [required] 3/ Output file name \n
       [required] 4/ An optional file containin blaclisted regions (tRNAs...) \n")
}

suppressMessages(library(data.table))

# Tests ----
# umi_count <- "db/counts/PROseq/HCFC1/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_counts.txt"
# annotation <- "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_promoters.rds"
# outputFile <- "db/count_tables/PROseq/HCFC1/promoter/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_promoter_counts.txt"

# Args ----
umi_count <- args[1]
annotation <- args[2]
outputFile <- args[3]
blacklist <- if(length(args)==4)
  args[4] else
    NULL

# Import annotation ----
annot <- readRDS(annotation)

# Import UMI ----
dat <- fread(umi_count)
dat[, c("seqnames", "start", "strand"):= tstrsplit(coor, ":", type.convert = T)]

# Remove reads overlapping blacklisted regions ----
if(!is.null(blacklist)) {
  # Import blacklisted regions
  blacklist <- readRDS(blacklist)
  blacklist <- data.table::copy(blacklist)
  # Compute overlap
  ov <- blacklist[dat, .N, on= c("seqnames", "start<=start", "end>=start", "strand"), .EACHI]$N
  # Print statistics
  stats <- data.table(
    N.removed.reads.clusters= sum(ov>0, na.rm = T),
    N.removed.umi.counts= sum(dat[ov>0, umi_counts], na.rm = T),
    perc.removed.umi.counts= sum(dat[ov>0, umi_counts], na.rm = T)/sum(dat$umi_counts, na.rm = T)*100
  )
  fwrite(stats,
         gsub(".txt$", "__blacklisted.txt", outputFile),
         sep= "\t",
         na= NA)
  # Remove overlaps
  dat <- dat[ov==0]
}

# Compute counts ----
annot$count <- dat[annot, sum(umi_counts, na.rm= T), on= c("seqnames", "start>=start", "start<=end", "strand"), .EACHI]$V1
fwrite(annot[, .(ID= cluster.id, count)],
       outputFile,
       sep= "\t",
       na= NA)
