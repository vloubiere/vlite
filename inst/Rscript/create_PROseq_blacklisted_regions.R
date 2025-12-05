setwd("/groups/stark/vloubiere/genomes/Mus_musculus/PROseq/")

# Import blaclisted regions UCSC
blacklist <- fread("mm10-blacklist.v2.bed.gz")
setnames(blacklist, c("seqnames", "start", "end", "name"))
blacklist[, start:= start+1]
blacklist <- rbind(
  data.table::copy(blacklist)[, strand:= "+"],
  data.table::copy(blacklist)[, strand:= "-"]
)

# Import tRNAs gtf file
tRNAs <- rtracklayer::import("/groups/stark/vloubiere/genomes/Mus_musculus/GENCODE/gencode.vM25.tRNAs.gtf.gz")
tRNAs <- as.data.table(tRNAs)

# Import genes gtf file
others <- rtracklayer::import("/groups/stark/vloubiere/genomes/Mus_musculus/GENCODE/gencode.vM25.annotation.gtf.gz")
others <- as.data.table(others)
other_mask_types <- c(
  "rRNA","Mt_rRNA",
  "Mt_tRNA",
  "snRNA","snoRNA","scaRNA",
  "miRNA",
  "scRNA","sRNA",
  "misc_RNA",
  "ribozyme"
)
others <- others[gene_type %in% other_mask_types]

# Combine
mask <- rbind(others[, .(seqnames, start, end, strand, cluster.id= gene_type)],
              tRNAs[, .(seqnames, start, end, strand, cluster.id= gene_type)],
              blacklist[, .(seqnames, start, end, strand, cluster.id= name)])
mask[, start:= as.integer(start)]
mask[, cluster.id:= paste0(cluster.id, "__", seqnames, ":", start, "-", end, ":", strand)]
mask[, width:= as.integer(end-start+1)]

# Order
setcolorder(mask, "width", after= "end")
setorderv(mask, c("seqnames", "start", "end"))

# SAVE
exportBed(mask, "mm10_blacklisted_regions.bed")
saveRDS(mask, "mm10_blacklisted_regions.rds")

# # Import ncRNA regions UCSC
# fa.seq <- seqinr::read.fasta("/groups/stark/vloubiere/genomes/Mus_musculus/ncRNAs/Mus_musculus.GRCm38.ncrna.fa")
# ncRNA <- data.table(
#   seq= unlist(tstrsplit(sapply(fa.seq, attr, "Annot"), " ", keep= 3)),
#   cluster.id= unlist(tstrsplit(sapply(fa.seq, attr, "Annot"), " ", keep= 5))
# )
# ncRNA <- ncRNA[grepl("^chromosome:GRCm38:", seq)]
# ncRNA[, cluster.id:= tstrsplit(cluster.id, ":", keep= 2)]
# ncRNA[, c("seqnames", "start", "end", "strand"):= tstrsplit(seq, ":", keep= 3:6, type.convert = T)]
# ncRNA <- ncRNA[!grepl("^CHR", seqnames)]
# ncRNA[, seqnames:= paste0("chr", seqnames)]
# ncRNA[, start:= start+1]
# ncRNA[, strand:= ifelse(strand==-1, "-", "+")]
# ncRNA[, width:= end-start+1]
# ncRNA$seq <- NULL
# ncRNA[, cluster.id:= paste0(cluster.id, "__", seqnames, ":", start, "-", end, ":", strand)]
# setcolorder(ncRNA, c("seqnames", "start", "end", "width", "strand", "cluster.id"))
