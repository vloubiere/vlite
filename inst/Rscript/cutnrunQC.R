#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4, 5)) {
  stop("Please specify:\n",
       "[required] 1/ Peaks file \n",
       "[required] 2/ Output pdf file \n",
       "[required] 3/ BSgenome name \n",
       "[required] 4/ Sample bigwig file \n",
       "[optional] 5/ Input bigwig file \n")
}

# Load packges ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(devtools, warn.conflicts = FALSE))
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Tests
# peaks <- "/groups/stark/nemcko/Hcfc1_paper/Hcfc1/db/peaks/CUTNRUN/Hcfc1.C_AID.Hcfc1_noTreatment_exp4_rep1_peaks.narrowPeak"
# peaks <- "/groups/stark/nemcko/Hcfc1_paper/Hcfc1/db/peaks/CUTNRUN/Hcfc1.C_AID.Hcfc1_IAA_6h_exp4_merged_peaks.narrowPeak"
# genome <- "mm10"
# bw <- "/groups/stark/nemcko/Hcfc1_paper/Hcfc1/db/bw/CUTNRUN/Hcfc1.C_AID.Hcfc1_noTreatment_exp4_rep1.bw"
# input <- "/groups/stark/nemcko/Hcfc1_paper/Hcfc1/db/bw/CUTNRUN/V5_parental_noTreatment_exp1_rep1.bw"

# Parse args
peaks <- args[1]
pdf.file <- args[2]
genome <- args[3]
bw <- args[4]
input <- if(length(args)==5)
  args[5] else
    NULL

# Import peaks ----
peaks <- importBed(peaks)

# Compute quantiles ----
peaks[, quant:= cut(signalValue,
                    quantile(signalValue,
                             seq(0, 1, length.out= 6)),
                    include.lowest = T)]
setorderv(peaks, "quant")

# Random control peaks ----
rdm <- vlite::randomRegionsBSgenome(genome = genome,
                                    widths = peaks[, end-start+1],
                                    restrict.seqnames = unique(peaks$seqnames),
                                    no.overlaps = peaks)

# Combine objects for average tracks ----
cmb <- rbind(rdm, peaks, fill= T)
by <- c(rep("rdm", nrow(rdm)), as.character(peaks$quant))
by <- factor(by, unique(by))

# Compute signal ----
pl <- if(is.null(input)) {
  list(rdm.sample= bwCoverage(rdm, bw),
       peaks.sample= bwCoverage(peaks, bw))
} else {
  list(rdm.input= bwCoverage(rdm, input),
       peaks.input= bwCoverage(peaks, input),
       rdm.sample= bwCoverage(rdm, bw),
       peaks.sample= bwCoverage(peaks, bw))
}

# Color legend ----
Cc <- colorRampPalette(c("pink", "red"))(5)

# Plot QC pdf ----
pdf(pdf.file, width = 11, height = 3.25)
layout(matrix(1:4, ncol= 4),
       widths = c(1,1.3,1,1.3))
vl_par(mai= c(.9, .5, .5, .5),
       omi= c(0, 0, .25, 1.5))
# Quantif signal
vl_boxplot(pl,
           compute.pval = if(is.null(input)) list(c(1,2)) else list(c(1,2), c(1,3), c(2,4), c(3,4)),
           col= rep(c("lightgrey", "pink"), each= length(pl)/2),
           notch= TRUE,
           ylab= "Signal",
           tilt.names = T)
mtext(text= gsub(".bw$", "", basename(bw)), outer = T)
# Sorted signal
vl_plot(sort(peaks$signalValue),
        col= Cc[peaks$quant],
        cex= .6,
        ylab= "signalValue")
vl_legend("topleft",
          col= Cc,
          legend= levels(peaks$quant),
          pch= 16)
# Signal value
vl_boxplot(signalValue~quant,
           peaks,
           ylab= "signalValue",
           tilt.names= TRUE,
           col= Cc)
# Average tracks
bwAverageTrack(bed = cmb,
               by= by,
               tracks = bw,
               col.palette = c("grey", Cc),
               names = "")
dev.off()
