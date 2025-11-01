setwd("/groups/stark/vloubiere/projects/PROseq_pipeline/")
require(BSgenome.Mmusculus.UCSC.mm10)

# This script was part of Vanja's pipeline. The only modification I did is line 14, to get more informative IDs ----

# Test ----
### one transcript per gene - from new mm10 PRO-seq annotation
load(file = "/groups/stark/haberle/general.data/gene.annotation/mm10/PROseq/One.per.gene.CAGE.defined.promoters.near.annotated.TSSs.with.associated.nearest.annotated.transcript.RData")  # nonoverlap.f5.tss.near.ann
a <- one.per.gene.f5.tss.near.ann

# sort by chromosome and position
a <- a[order(a$chr, a$TSS),]
a$diff_analysis_id <- paste(a$chr, a$TSS, a$strand, sep = "_")

# >>>>>>> Modified by me
a$cluster.id <- paste0(a$gene_id, "__", a$gene_name, "__", a$diff_analysis_id)
# <<<<<<<

# define promoter and gene body regions (relative to the CAGE TSS associated with each transcript)
flank.up <- 0
flank.down <- 150
tss.gr <- GRanges(seqnames = a$chr,
                  ranges = IRanges(start = a$TSS, end = a$TSS),
                  strand = a$strand,
                  cluster.id = a$cluster.id,
                  seqlengths = seqlengths(Mmusculus))
promoter.gr <- promoters(tss.gr,
                         upstream=flank.up,
                         downstream=flank.down)
genebody.gr <- GRanges(seqnames = a$chr,
                       ranges = IRanges(start = ifelse(a$strand == "+", a$TSS+flank.down, a$transcript_start),
                                        end = ifelse(a$strand == "+", a$transcript_end, a$TSS-flank.down)),
                       strand = a$strand,
                       cluster.id = a$cluster.id)
transcript.gr <- GRanges(seqnames = a$chr,
                         ranges = IRanges(start = ifelse(a$strand == "+", a$TSS, a$transcript_start),
                                          end = ifelse(a$strand == "+", a$transcript_end, a$TSS)),
                         strand = a$strand,
                         cluster.id = a$cluster.id)

# Save files ----
saveRDS(as.data.table(promoter.gr), "db/annotations/mm10_promoters.rds")
saveRDS(as.data.table(genebody.gr), "db/annotations/mm10_genebody.rds")
saveRDS(as.data.table(transcript.gr), "db/annotations/mm10_transcript.rds")
