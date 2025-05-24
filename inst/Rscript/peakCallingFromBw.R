#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 10) {  # Fix the condition to match required arguments
  stop("Please specify:\n",
       "[required] 1/ Experiment bw file. \n",
       "[required] 2/ Optional input bw file. If left empty (''), local background only will be used. \n",
       "[required] 3/ Optional path to a bed file, specifying a subset of regions for which peaks should be called. \n",
       "[required] 4/ Output .narrowPeak file path. \n",
       "[required] 5/ Output .txt statistics file. \n",
       "[required] 6/ z-score cutoff to identify putative peaks. Default= 1.64. \n",
       "[required] 7/ Local background enrichment cutoff. Default= 1.5. \n",
       "[required] 8/ Local background FDR cutoff. Default= 0.05. \n",
       "[required] 9/ Input enrichment cutoff. Default= 1.5. \n",
       "[required] 10/ Input FDR cutoff. Default= 0.05. \n")
}

# Load packges ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomeInfoDb, warn.conflicts = FALSE))
devtools::load_all("/groups/stark/vloubiere/vlite/")

# # Tests ----
# track <- "/groups/stark/vloubiere/projects/Facilitators/db/bw/merge/screen_1500bp_K562_merged_mm9.bw"
# input <- "/groups/stark/vloubiere/projects/Facilitators/db/bw/merge/input_1500bp_none_merged_mm9.bw"
# bed.subset <- readRDS("/groups/stark/vloubiere/projects/Facilitators/Rdata/regions_of_interest.rds")[type=="locus"]
# zscore.cutoff <- 1.64
# local.enr.cutoff <- 1.5
# local.fdr.cutoff <- 0.05
# input.enr.cutoff <- 2
# input.fdr.cutoff <- 0.05

# Parse arguments ----
track <- args[1]
input <- if(args[2] == "NULL") NULL else args[2] # ifelse does not work here!
bed.subset <- if(args[3] == "NULL") NULL else args[3] # ifelse does not work here!
peak.file <- args[4]
stat.file <- args[5]
zscore.cutoff <- as.numeric(args[6])
local.enr.cutoff <- as.numeric(args[7])
local.fdr.cutoff <- as.numeric(args[8])
input.enr.cutoff <- as.numeric(args[9])
input.fdr.cutoff <- as.numeric(args[10])

# Checks ----
if(!is.null(bed.subset)) {
  stopifnot(grepl(".bed$", bed.subset))
  stopifnot(file.exists(bed.subset))
}
stopifnot(grepl(".narrowPeak$", peak.file))
stopifnot(grepl(".txt$", stat.file))

# Retrieve provided genome/regions coordinates ----
regions <- if(is.null(bed.subset)) vlite::bwGetSeqlengths(track) else vlite::importBed(bed.subset)
regions <- regions[, .(seqnames, start, end)]
regions[, region.idx:= .I]

# Save printed messages in stat.file ----
sink(stat.file, split = TRUE)

# Go for each region in genome (chromosomes or custom) ----
peaks <- regions[, {
  print(paste0("Starting with -> ", seqnames, ":", start, "-", end))
  t0 <- Sys.time()

  # Bin the whole region (step= 1) and compute signal ----
  bins <- vlite::binBed(.SD, bins.width = 100, steps.width = 1)
  bins[, signal:= vlite::bwCoverage(.SD, track)]

  # Subset non-overlapping bins (step= 100) to infer pseudocount, mean and sd ----
  non.overlapping <- bins[start %% 100 == 0]
  pseudo <- if(any(bins$signal==0))
    quantile(non.overlapping$signal[non.overlapping$signal > 0], 0.01) else
      0
  mu <- mean(log2(non.overlapping$signal+pseudo))
  sigma <- sd(log2(non.overlapping$signal+pseudo))

  # Identify candidate peaks using z-score ----
  bins[, zscore := (log2(signal+pseudo) - mu) / sigma]
  cand <- bins[signal>0 & zscore>zscore.cutoff]
  cand <- vlite::collapseBed(cand, min.gapwidth = 101)
  cand[, name:= paste0("peak_", .I, "_region_", region.idx)]

  # Compute signal in smaller bins (width= 50) with larger step (step= 25) ----
  large.bins <- vlite::binBed(.SD, bins.width = 50, steps.width = 25)
  large.bins[, signal:= vlite::bwCoverage(.SD, track)]
  if(!is.null(input)) large.bins[, input:= bwCoverage(.SD, input)]
  # Remove empty bins
  large.bins <- large.bins[signal>0]
  if(!is.null(input)) large.bins <- large.bins[input>0]
  # Normalize the signal
  large.bins[, signal:= signal/sum(signal)*1e6]
  if(!is.null(input)) large.bins[, input:= input/sum(input)*1e6]

  # Compute overlaps between large.bins and candidate peaks ----
  ov <- vlite::overlapBed(cand, large.bins, minoverlap = 40)
  # Remove candidate peaks with not overlaps
  cand <- cand[unique(ov$idx.a)]
  # Retrieve signal
  cand$signal <- ov[, large.bins[idx.b, .(.(signal))], keyby= idx.a]$V1
  if(!is.null(input))
    cand$input <- ov[, large.bins[idx.b, .(.(input))], keyby= idx.a]$V1

  # Select closest bins to use as local background ----
  bg.bins <- vlite::intersectBed(large.bins, cand, maxgap = 500, invert = T)
  closest <- vlite::closestBed(
    a = cand,
    b = bg.bins,
    k = max(lengths(cand$signal)) # Max number of bins within peak
  )

  # Sample closest local background bins to match the number of large.bins per peak ----
  closest[, Nbins:= lengths(cand$signal)[idx.a]]
  closest <- closest[order(abs(dist))]
  closest <- closest[, .SD[1:Nbins], .(idx.a, Nbins)]
  # Retrieve local bg signal from k closest non-overlapping bins
  cand$local.sig <- closest[, .(.(bg.bins$signal[idx.b])), keyby= idx.a]$V1

  # Return----
  t1 <- Sys.time()
  print(paste0("Done -> ", round(t1-t0, 2), "s"))
  cand
}, region.idx]
print(paste(nrow(peaks), "candidate peaks found."))

# Compute enrichment compared to local background ----
peaks[, c("signalValue", "pValue"):= {
  .(mean(signal[[1]])/mean(local.sig[[1]]),
    wilcox.test(signal[[1]],
                local.sig[[1]],
                alternative= "greater")$p.value)
}, name]
peaks[, qValue:= p.adjust(pValue, method= "fdr")]
peaks[, hit:= signalValue>=local.enr.cutoff & qValue<=local.fdr.cutoff]
peaks <- peaks[(hit)]
print(paste(nrow(peaks), "candidate peaks were enriched vs. local background."))

# Compute enrichment compared to input ----
if(!is.null(input) && nrow(peaks)) {
  peaks[, c("signalValue", "pValue"):= {
    .(mean(signal[[1]])/mean(input[[1]]),
      wilcox.test(signal[[1]],
                  input[[1]],
                  alternative= "greater")$p.value)
  }, name]
  peaks[, qValue:= p.adjust(pValue, method= "fdr")]
  peaks[, hit:= signalValue>=input.enr.cutoff & qValue<=input.fdr.cutoff]
  peaks <- peaks[(hit)]
  print(paste(nrow(peaks), "candidate peaks were enriched vs. input."))
}

if(nrow(peaks)) {
  # Final peaks name ----
  peaks[, name:= paste0("peak_", .I)]

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
} else {
  peaks[, score:= numeric()]
  peaks[, peak:= integer()]
}

# Save peaks and statistics file ----
exportBed(peaks, peak.file)
sink()
