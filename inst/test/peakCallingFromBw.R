setwd("/groups/stark/vloubiere/vlite/")
devtools::load_all("./")

old <- importBed("/groups/stark/vloubiere/projects/epigenetic_cancer/db/peaks/cutnrun/H3K27Ac_PH18_confident_peaks.narrowPeak")
old <- old[seqnames=="chr2L"]

track <- "/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/cutnrun/H3K27Ac_PH18_merge.bw"
# chr <- importBed("chr2L:4000000-5000000")
chr <- importBed("chr2L:1-23513712")
# Compute bins
bins <- binBed(chr, bins.width = 100, steps.width = 1)
bins <- bins[, .(seqnames, start, end)]

# Compute signal
bins[, signal:= bwCoverage(.SD, track)]

# Remove bins with no signal
bins <- bins[signal>0]

# Log2 transform
bins[, signal:= log2(signal)]

# Gaussian blue?
# bins[, signal:= vlite::gaussianBlur(signal)]
# plot(density(bins$signal))

# Compute z-score at 3 different scales?
bins[, chr:= scale(signal)>1.64]

# Add 50000 bins (1/2 the large window, see below) to avoid NAs at start/end
start <- bins[, .(start= NA, end= NA, signal= mean(signal[1:5000]), chr= F), seqnames]
end <- bins[, .(start= NA, end= NA, signal= mean(signal[(.N-4999):.N]), chr= F), seqnames]
bins <- rbind(start[rep(1, 50000)], bins, end[rep(1, 50000)])

# Are there at least two significant z-scores within 100kb windows?
dists <- c(large.window= 1e5, short.window= 1e4)
for(i in seq(dists)) {
  dist <- dists[i]
  line.s <- seq(1, dist+1, length.out= 11)[-11]
  for(j in seq(line.s)) {
    bins[line.s[j]:.N, line.idx:= round(.I/dist)]
    bins[line.s[j]:.N, paste0("p.", j):= scale(signal)>1.64, line.idx]
  }
  bins[, names(dist):= rowSums(.SD)>1, .SDcols= patterns("^p.")]
  bins <- bins[, !grep("^p.", names(bins), value = T), with= F]
}

# Collapse putative peaks
peaks <- bins[(chr & large.window & short.window)]
peaks[, peak.idx:= collapseBed(.SD, return.idx.only = T, max.gap = 50)]
peaks[, center:= which.max(signal)-1, peak.idx]
peaks <- peaks[, .(start= start[1], end= end[.N]), .(seqnames, peak.idx, center)]

# Bin putative peaks in 25bp bins
peaks.bins <- binBed(peaks, bins.width = 25, steps.width = 25, bins.width.min = 20)
peaks.bins[, signal:= bwCoverage(.SD, track)]
peaks$signal <- peaks.bins[, .(.(signal)), peak.idx]$V1

# Subtract putative peaks from full chromosome
noPeaks <- subtractBed(chr, peaks)
noPeaks <- binBed(noPeaks, bins.width = 25, steps.width = 25, bins.width.min = 25)
noPeaks[, signal:= bwCoverage(.SD, track)]
noPeaks <- intersectBed(noPeaks,
                        resizeBed(peaks,
                                  center = "center",
                                  upstream= 10000,
                                  downstream= 10000,
                                  genome = "dm6"))

# For each putative peak, identify the 20 closest bins
cl <- closestBed(peaks, noPeaks, n=20, min.dist = 1)

# Retrieve signal at control bins
cl[, ctl.signal:= noPeaks$signal[idx.b]]
peaks$ctl.signal <- cl[, .(.(ctl.signal)), idx.a]$V1

# Compute p.value
peaks[, pval:= wilcox.test(signal[[1]], ctl.signal[[1]], alternative= "greater")$p.value, peak.idx]
peaks[, padj:= p.adjust(pval, method = "fdr")]
plot(sort(-log10(peaks$padj)))

# Compute Enrichment
peaks[, enr:= mean(signal[[1]])/mean(ctl.signal[[1]]), peak.idx]

# Select significant peaks
plot(peaks$enr, -log10(peaks$padj))
sel <- peaks[padj<0.05 & enr>1.5]

# Save bed file
