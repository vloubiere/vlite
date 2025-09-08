setwd("/groups/stark/vloubiere/projects/scRNASeq_ED/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import motifs ----
dat <- vlite::importJASPAR("/groups/stark/vloubiere/motifs_db/hand_curated_Dmel_motifs_SCENIC_lite.txt")

# Import FPKM data ----
FPKM <- data.table(file= list.files("db/FPKM/RNA/", full.names = T))
FPKM[, cdition:= tstrsplit(basename(file), "__", keep= 1)]
FPKM <- FPKM[, fread(file), cdition]

# Compute TF metadata ----
meta <- data.table(motif= names(dat$pwms_log_odds))
meta[, gene_symbol:= tstrsplit(motif, "__", keep= 1)]

# Add composite motifs ----
compMot <- grep("::", meta$gene_symbol, value= T)
for(comp in compMot) {
  .c <- unlist(tstrsplit(comp, "::"))
  .c <- FPKM[gene_symbol %in% .c, lapply(.SD, min), cdition, .SDcols= c("mean_FPKM", "log2_FPKM")]
  .c[, gene_symbol:= comp]
  FPKM <- rbind(FPKM, .c, fill= T)
}
meta <- merge(meta, FPKM, by= "gene_symbol")
setcolorder(meta, "motif")

# Save final object ----
final <- c(list(meta= meta), dat)
saveRDS(final,
        "/groups/stark/vloubiere/motifs_db/Dmel_motifs_SCENIC_lite_db.rds")
