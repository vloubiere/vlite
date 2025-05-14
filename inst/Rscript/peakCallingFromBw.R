#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly = TRUE)
if(!length(args) == 11) {  # Fix the condition to match required arguments
  stop("Please specify:\n",
       "[required] 1/ Experiment bw file. \n",
       "[required] 2/ Optional input bw file. If left empty (''), local bacground only will be used. \n",
       "[required] 3/ Path to a bed file containing the regions from which peaks should be called. \n",
       "[required] 4/ A BSgenome name for which peaks should be called. If specified, overrides bed argument (3) \n",
       "[required] 5/ Output .narrowPeak file path. \n",
       "[required] 6/ Output .txt statistics file. \n",
       "[required] 7/ z-score cutoff to identify putative peaks. Default= 1.64. \n",
       "[required] 8/ Local background enrichment cutoff. Default= 1.5. \n",
       "[required] 9/ Local background FDR cutoff. Default= 0.05. \n",
       "[required] 10/ Input enrichment cutoff. Default= 1.5. \n",
       "[required] 11/ Input FDR cutoff. Default= 0.05. \n")
}

# Load packges ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomeInfoDb, warn.conflicts = FALSE))
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Parse arguments ----
track <- args[1]
input <- if(args[2] == "NULL") NULL else args[2] # ifelse does not work here!
bed <- if(args[3] == "NULL") NULL else args[3] # ifelse does not work here!
genome <- if(args[4] == "NULL") NULL else args[4] # ifelse does not work here!
peak.file <- args[5]
stat.file <- args[6]
zscore.cutoff <- as.numeric(args[7])
local.enr.cutoff <- as.numeric(args[8])
local.fdr.cutoff <- as.numeric(args[9])
input.enr.cutoff <- as.numeric(args[10])
input.fdr.cutoff <- as.numeric(args[11])

# Retrieve genome regions ----
if(!is.null(genome)) {
  bed <- tempfile(fileext = ".bed")
  vlite::exportBed(bed = vlite::getBSgenomeSize(genome),
                   file = regions)

}
if(is.null(bed))
  stop("Bed file could not be found")

# Checks ----
if(!grepl(".bed$", bed)) stop("bed should be a path to a unique bed file.")
if(!file.exists(bed)) stop("bed file could not be found.")
if(!grepl(".narrowPeak$", peak.file)) stop("peak.file should end up with .narrowPeak extension")
if(!grepl(".txt$", stat.file)) stop("stat.file should end up with .txt extension")

# Save logs ----
sink(stat.file, split = TRUE)

# Retrieve provided genome/regions coordinates ----
chr <- vlite::bwGetSeqnames(track)
genome <- vlite::importBed(bed)
genome <- genome[seqnames %in% chr, .(seqnames, start, end)]
genome[, region.idx:= .I]

# # Tests ----
# track <- "/groups/stark/vloubiere/projects/Facilitators/db/bw/merge/screen_1500bp_K562_merged_mm9.bw"
# input <- "/groups/stark/vloubiere/projects/Facilitators/db/bw/merge/input_1500bp_none_merged_mm9.bw"
# genome <- readRDS("/groups/stark/vloubiere/projects/Facilitators/Rdata/regions_of_interest.rds")[type=="locus"]
# zscore.cutoff <- 1.64
# local.enr.cutoff <- 1.5
# local.fdr.cutoff <- 0.05
# input.enr.cutoff <- 2
# input.fdr.cutoff <- 0.05

# Go for each region in genome (chromosome or custom) ----
peaks <- genome[, {
  print(paste0("Starting with -> ", seqnames, ":", start, "-", end))
  t0 <- Sys.time()
  # Identify putative peaks using zscore ----
  # Bin Region
  bins <- binBed(.SD, bins.width = 100, steps.width = 1)
  # Compute signal
  bins[, signal:= bwCoverage(.SD, track)]
  # Subset non-overlapping bins
  non.overlapping <- bins[start %% 100 == 0]
  # Compute pseudocount
  pseudo <- quantile(non.overlapping$signal[non.overlapping$signal > 0], 0.01)
  # Compute mean and sd
  mu <- mean(log2(non.overlapping$signal+pseudo))
  sigma <- sd(log2(non.overlapping$signal+pseudo))
  # Compute z-score
  bins[, zscore := (log2(signal+pseudo) - mu) / sigma]
  # Retrieve putative peaks
  cand <- bins[signal>0 & zscore>zscore.cutoff]
  cand <- collapseBed(cand, min.gapwidth = 101)
  cand[, name:= paste0("peak_", collapseBed(.SD, return.idx.only = T))]

  # Test peaks enrichment using smaller bins with larger step ----
  # Compute bins with large step
  large.bins <- binBed(.SD, bins.width = 50, steps.width = 25)
  # Compute signal
  large.bins[, signal:= bwCoverage(.SD, track)]
  if(!is.null(input)) large.bins[, input:= bwCoverage(.SD, input)]
  # Remove empty bins
  large.bins <- large.bins[signal>0]
  if(!is.null(input)) large.bins <- large.bins[input>0]
  # Normalize the signal
  large.bins[, signal:= signal/sum(signal)*1e6]
  if(!is.null(input)) large.bins[, input:= input/sum(input)*1e6]

  # Retrieve signal for cancidate peaks ----
  # Compute overlap large.bins / candidate peaks
  ov <- overlapBed(cand, large.bins, minoverlap = 40)
  # Remove candidate peaks with not overlaps
  cand <- cand[unique(ov$idx.a)]
  # Add signal to candidate peaks
  cand$signal <- ov[, large.bins[idx.b, .(.(signal))], keyby= idx.a]$V1
  if(!is.null(input))
    cand$input <- ov[, large.bins[idx.b, .(.(input))], keyby= idx.a]$V1
  # Select closest bins to use as local background
  bg.bins <- intersectBed(large.bins, cand, maxgap = 500, invert = T)
  closest <- closestBed(
    a = cand,
    b = bg.bins,
    k = max(lengths(cand$signal)) # Max number of bins within peak
  )
  # Subset local background bins to match the number of bins within peaks
  closest[, Nbins:= lengths(cand$signal)[idx.a]]
  closest <- closest[order(abs(dist))]
  closest <- closest[, .SD[1:Nbins], .(idx.a, Nbins)]
  # Retrieve local bg signal from k closest non-overlapping bins
  cand$local.sig <- closest[, .(.(bg.bins$signal[idx.b])), keyby= idx.a]$V1
  # Return
  t1 <- Sys.time()
  print(paste0("Done -> ", round(t1-t0, 2), "s"))
  cand
}, region.idx]
print(paste(nrow(peaks), "candidate peaks found."))

# Compute enrichment compared to local background
peaks[, c("signalValue", "pValue"):= {
  .(mean(signal[[1]])/mean(local.sig[[1]]),
    wilcox.test(signal[[1]],
                local.sig[[1]],
                alternative= "greater")$p.value)
}, name]
peaks[, qValue:= p.adjust(pValue, method= "fdr")]
peaks[, hit:= signalValue>local.enr.cutoff & qValue<local.fdr.cutoff]
peaks <- peaks[(hit)]
print(paste(nrow(peaks), "candidate peaks were enriched vs. local background."))

# Compute enrichment compared to input
if(!is.null(input)) {
  peaks[, c("signalValue", "pValue"):= {
    .(mean(signal[[1]])/mean(input[[1]]),
      wilcox.test(signal[[1]],
                  input[[1]],
                  alternative= "greater")$p.value)
  }, name]
  peaks[, qValue:= p.adjust(pValue, method= "fdr")]
  peaks[, hit:= signalValue>input.enr.cutoff & qValue<input.fdr.cutoff]
  peaks <- peaks[(hit)]
  print(paste(nrow(peaks), "candidate peaks were enriched vs. input."))
}

# Select peaks and compute score ----
peaks[, score:= floor(-log10(qValue) / max(-log10(qValue)) * 1000)]
peaks[score>1000, score:= 1000]

# Find summit ----
t0 <- Sys.time()
summit <- binBed(peaks[, .(seqnames, start, end)], bins.width = 1, steps.width = 1)
summit[, signal:= bwCoverage(.SD, track)]
peaks$peak <- summit[, .(which.max(signal)-1), keyby= line.idx]$V1
t1 <- Sys.time()
print(paste0("Finding summits -> ", round(t1-t0, 2), "s"))

# Save peaks and statistics file ----
exportBed(peaks, peak.file)
sink()
