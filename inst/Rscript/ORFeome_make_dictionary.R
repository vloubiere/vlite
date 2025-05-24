#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# README
# Confident BC/ORFs pairs are defined as follows:
# 1/ BCs that have the expected lenght and sequence pattern are selected
# 2/ Uniquely align ORF reads with a mapping quality >= 30 are selected
# 3/ BCs/ORFs reads are paired based on read ID, and only BCs associated to a unique ORF are retained
# 4/ Finally, only BC/ORF reads supported by at least x independent reads are retained (see arg 3 below)

# Check arguments ----
if (length(args) != 4) {
  stop("Please specify:\n
       [required] 1/ A comma-separated list of fq.gz files containing the BC sequences \n
       [required] 2/ A bam file containing ORF alignments for all the reads present in the BC file (arg[1]) \n
       [required] 3/ Minimum number of supporting reads \n
       [required] 4/ Output prefix for output files")
}
require(Rsamtools)
require(rtracklayer)
require(GenomicRanges)
require(data.table)
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Variables ----
trim_BC <- strsplit(args[1], ",")[[1]]
bam <- args[2]
minNumberReads <- as.numeric(args[3])
output_prefix <- args[4]
stopifnot(all(grepl(".fq.gz$", trim_BC)))
stopifnot(grepl(".bam$", bam))
stopifnot(minNumberReads %% 1 == 0)

# Tests ----
# meta <- readRDS("Rdata/processed_dictionary_metadata.rds")[sampleID=="lib200S16_dictionary"]
# bam <- unique(meta$bam)
# trim_BC <- unique(meta$trim_BC)
# minNumberReads <- 3
# output_prefix <- "lib200S16_dictionary"

# Import BC files ----
BC <- lapply(trim_BC, vlite::importFq)
BC <- rbindlist(BC)
setnames(BC, "Sequence", "BC")
BC <- BC[!is.na(BC)]
# Quality check (length between 25 and 30 + matches expected pattern)
BC[, BC:= as.character(reverseComplement(DNAStringSet(BC)))]
BC[, qual:= grepl("([GC][AT]){4}[GCAT]{5}([GC][AT]){5}[GCAT]{2}", BC) & between(nchar(BC), 25, 30) , BC]
# Filter BC quality
BCq <- BC[(qual)]

# Import ORFs bam file ----
ORF <- vlite::importBam(bam,
                        sel= c("qname", "rname", "mapq", "strand"),
                        new.col.names = c("readID", "seqnames", "mapq", "strand"))
ORF[, seqnames:= as.character(seqnames)]
# Select best alignment per readID
setorderv(ORF, "mapq", -1, na.last = T)
ORF <- ORF[, .SD[1], readID]
# Filter mapq quality
ORFq <- ORF[mapq>=30]

# Connect BCs to the ORFs ----
pairs <- merge(BCq, ORFq, by= "readID")
pairs[, total:= .N, BC]
pairs[, count:= .N, .(seqnames, BC)]

# Collapse pairs ----
coll <- unique(pairs[, .(BC, seqnames, total, count)])

# Filter usable BCs, with >= minNumberReads reads that are all assigned to the same ORF ----
final <- coll[count==total & count>=minNumberReads]

# Calculate read-level statistics ----
read_stats <- c(
  "Valid BC+ORF" = nrow(pairs),
  "Valid ORF"    = nrow(ORFq) - nrow(pairs),
  "Valid BC"     = nrow(BCq) - nrow(pairs)
)
read_stats["None"] <- nrow(ORF) - sum(read_stats)

# Calculate barcode-level statistics (using only valid reads) ----
bc_stats_values <- c(
  nrow(final),
  nrow(coll[count >= minNumberReads]) - nrow(final),
  nrow(coll[count == total]) - nrow(final)
)
bc_stats_names <- c(
  paste0("Uniq.ORF + counts>=", minNumberReads),
  paste0("Counts>=", minNumberReads),
  "Uniq.ORF"
)
names(bc_stats_values) <- bc_stats_names
bc_stats <- bc_stats_values
bc_stats["None"] <- nrow(coll) - sum(bc_stats)

# Compute barcodes per ORF ----
hdr <- Rsamtools::scanBamHeader(bam)
all_seqnames <- data.table(seqnames = names(hdr[[1]]$targets))
pairs_stats <- final[all_seqnames, .N, .EACHI, on = "seqnames"]

# Prepare summary statistics ----
stats <- c(
  sprintf("Total reads: %d", nrow(ORF)),
  "Of these:",
  paste0("\t", names(read_stats), ": ", read_stats),
  sprintf("Usable reads (Valid BC+ORF) contained: %d in total", nrow(coll)),
  "Of these:",
  paste0("\t", names(bc_stats), ": ", bc_stats),
  sprintf(
    "Out of %d ORFs in total, %d had at least 5 BCs",
    nrow(pairs_stats),
    sum(pairs_stats$N >= 5)
  )
)

# Save statistics ----
writeLines(stats, paste0(output_prefix, "_stats.txt"))

# Plot statistics ----
pdf(paste0(output_prefix, "_stats.pdf"), 9, 3)
vl_par(mfrow= c(1,3))
# Read counts 
pie(read_stats,
    labels = paste0(names(read_stats), "\n(n=", formatC(read_stats, big.mark = ","), ")"),
    main = paste("Read counts", basename(output_prefix)))
text(x= 0.5,
     y= par("usr")[3],
     paste0("Continuing with ", formatC(read_stats[1], big.mark = ","), " usable reads"),
     xpd= NA,
     pos= 1,
     offset= 1)
# ORF counts 
pie(bc_stats,
    labels = paste0(names(bc_stats), "\n(n=", formatC(bc_stats, big.mark = ","), ")"),
    main= paste("BC counts", basename(output_prefix)))
text(x= 0.5,
     y= par("usr")[3],
     paste0("BCs with uniq.ORF & counts>=5 \nrepresent ",
            round(sum(final$count)/read_stats[1]*100), "% of usable reads"),
     xpd= NA,
     pos= 1,
     offset= 1)
# Retrieve all sequence in bam index
clipped.N <- ifelse(pairs_stats$N>500, 500, pairs_stats$N)
hist(clipped.N,
     breaks = seq(0, max(clipped.N)),
     border= NA,
     col= "black",
     xlab= "Number of unambiguous BCs",
     main= paste("Assigned BCs", basename(output_prefix)))
vl_legend("topright",
          legend= paste0(sum(pairs_stats$N>5), "/", nrow(pairs_stats),
                         " ORFs\nwith >= 5 BCs"))
dev.off()

# Save
saveRDS(final,
        paste0(output_prefix, ".rds"))