# Download Dev/Hk enhancers from pe-STARR-Seq paper
tmp <- tempfile(pattern = ".xlsx")
download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-52921-2/MediaObjects/41467_2024_52921_MOESM4_ESM.xlsx",
             destfile = tmp)

# Retrieve enhancers
enh <- readxl::read_xlsx(tmp, sheet = 2, skip = 1)
enh <- as.data.table(enh)
enh <- enh[group=="dev" & detail %in% c("medium", "strong") | group=="hk"]

# Negative controls
set.seed(1)
ctl <- controlRegionsBSgenome(bed = enh,
                              genome = "dm3",
                              no.overlap = TRUE)

# Select JASPAR motifs
load("/groups/stark/vloubiere/motifs_db/vl_Dmel_motifs_DB_full.RData")
sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds]
pwms <- vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_perc]
mot.cluster <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_cluster]

# Compute enrichment at developmental enhancers
counts <- vl_motifCounts(bed= enh[group=="dev"],
                         genome= "dm3",
                         pwm_log_odds= sel)
ctl.counts <- vl_motifCounts(bed= ctl,
                             genome= "dm3",
                             pwm_log_odds= sel)

# Enrichment
enr <- vl_motifEnrich(counts,
                      control.counts = ctl.counts,
                      names = mot.cluster)

# Collapse per motif cluster
coll <- enr[motif %in% enr[, motif[which.min(padj)], name]$V1]

# Plot
vl_par(mai= c(.9, 3, .9, 1.2))
pl <- plot(coll)
addMotifs(plot.DT = pl,
          pwms = pwms)

# Compare dev and Hk enhancers
counts1 <- vl_motifCounts(bed= enh,
                          genome= "dm3",
                          pwm_log_odds= sel)

# Enrichment per cluster
enr1 <- vl_motifEnrich(split(counts1, enh$group),
                       control.counts = ctl.counts,
                       names = mot.cluster)

# Collapse per motif cluster and plot
coll1 <- enr1[motif %in% enr1[, motif[which.min(padj)], name]$V1]

# Plot
vl_par(mai= c(.9, 3, .9, 1.2))
pl <- plot(coll1)
addMotifs(plot.DT = pl,
          pwms= pwms)

