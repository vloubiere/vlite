# Import tRNAs gtf file
tRNAs <- rtracklayer::import("/groups/stark/vloubiere/genomes/Mus_musculus/GENCODE/gencode.vM25.tRNAs.gtf.gz")
tRNAs <- as.data.table(tRNAs)

# Import genes gtf file
genes <- rtracklayer::import("/groups/stark/vloubiere/genomes/Mus_musculus/GENCODE/gencode.vM25.annotation.gtf.gz")
genes <- as.data.table(genes)
genes_mask_types <- c(
  "rRNA","Mt_rRNA",
  "Mt_tRNA",
  "snRNA","snoRNA","scaRNA",
  "miRNA",
  "scRNA","sRNA",
  "misc_RNA",
  "ribozyme"
)
genes <- genes[gene_type %in% genes_mask_types]

# Combine
mask <- rbind(genes[, .(seqnames, start, end, gene_type, strand)], tRNAs[, .(seqnames, start, end, strand, gene_type)])
mask[, idx:= collapseBed(mask, return.idx.only = T, ignore.strand = FALSE)]
mask[, gene_type:= paste0(sort(unique(gene_type)), collapse= ","), idx]
mask[, start:= min(start), idx]
mask[, end:= max(end), idx]
mask <- unique(mask)
mask$idx <- NULL
mask[, width:= end-start+1]
mask[, cluster.id:= paste0(gene_type, "__", seqnames, ":", start, "-", end, ":", strand)]
mask$gene_type <- NULL
setcolorder(mask, "width", after= "end")

saveRDS(mask, "/groups/stark/vloubiere/genomes/Mus_musculus/PROseq/mm10_blacklisted_regions.rds")
