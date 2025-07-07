#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Test if there are 12 args: if not, return an error
if (length(args)!=12) {
  stop("Please specify:\n
       [required] 1/ A comma-separated list of count (ref genome)\n
       [required] 2/ A comma-separated list of read statistics (reference genome) \n
       [required] 3/ A comma-separated list of spike-in statistics \n
       [required] 4/ A comma-separated list of sample names \n
       [required] 5/ A comma-separated list of conditions \n
       [required] 6/ A comma-separated list of controls \n
       [required] 7/ dds output folder \n
       [required] 8/ FC tables output folder \n
       [required] 9/ PDF output folder \n
       [required] 10/ Experiment \n
       [required] 11/ feature \n
       [required] 12/ Normalization method. Possible values are 'default' (DESeq2 default), 'libSize', 'spikeIn', 'combined' (Vanja's method, for whichchanging the control sample will change the outcome)\n")
}

# Load libraries
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(DESeq2, warn.conflicts = FALSE))

# Tests ----
# counts <- c("/groups/stark/vloubiere/projects/vl_pipelines/db/count_tables/PROseq/HCFC1/transcript/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_transcript_counts.txt",
#             "/groups/stark/vloubiere/projects/vl_pipelines/db/count_tables/PROseq/HCFC1/transcript/AID-Hcfc1-cl17_0hrIAA_HCFC1_rep1_mm10_transcript_counts.txt",
#             "/groups/stark/vloubiere/projects/vl_pipelines/db/count_tables/PROseq/HCFC1/transcript/AID-Hcfc1-cl4_3hrIAA_HCFC1_rep1_mm10_transcript_counts.txt",
#             "/groups/stark/vloubiere/projects/vl_pipelines/db/count_tables/PROseq/HCFC1/transcript/AID-Hcfc1-cl17_3hrIAA_HCFC1_rep1_mm10_transcript_counts.txt")
# refStats <- c("/groups/stark/vloubiere/projects/vl_pipelines/db/alignment_stats/PROseq/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_statistics.txt",
#               "/groups/stark/vloubiere/projects/vl_pipelines/db/alignment_stats/PROseq/AID-Hcfc1-cl17_0hrIAA_HCFC1_rep1_mm10_statistics.txt",
#               "/groups/stark/vloubiere/projects/vl_pipelines/db/alignment_stats/PROseq/AID-Hcfc1-cl4_3hrIAA_HCFC1_rep1_mm10_statistics.txt",
#               "/groups/stark/vloubiere/projects/vl_pipelines/db/alignment_stats/PROseq/AID-Hcfc1-cl17_3hrIAA_HCFC1_rep1_mm10_statistics.txt")
# spikeStats <- c("/groups/stark/vloubiere/projects/vl_pipelines/db/alignment_stats/PROseq/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_dm3_spikeIn_statistics.txt",
#                 "/groups/stark/vloubiere/projects/vl_pipelines/db/alignment_stats/PROseq/AID-Hcfc1-cl17_0hrIAA_HCFC1_rep1_dm3_spikeIn_statistics.txt",
#                 "/groups/stark/vloubiere/projects/vl_pipelines/db/alignment_stats/PROseq/AID-Hcfc1-cl4_3hrIAA_HCFC1_rep1_dm3_spikeIn_statistics.txt",
#                 "/groups/stark/vloubiere/projects/vl_pipelines/db/alignment_stats/PROseq/AID-Hcfc1-cl17_3hrIAA_HCFC1_rep1_dm3_spikeIn_statistics.txt")
# names <- c("control_rep1","control_rep2","IAA_3h_rep1","IAA_3h_rep2")
# conditions <- c("control","control","IAA_3h","IAA_3h")
# controls <- "control"
# dds_output_folder <- "db/dds/PROseq/"
# FC_output_folder <- "db/FC_tables/PROseq/"
# PDF_output_folder <- "pdf/PROseq/"
# experiment <- "HCFC1"
# feature <- "transcript"
# norm <- "spikeIn"

# Parse variables ----
counts <- unlist(tstrsplit(args[1], ","))
refStats <- unlist(tstrsplit(args[2], ","))
spikeStats <- unlist(tstrsplit(args[3], ","))
names <- unlist(tstrsplit(args[4], ","))
names <- gsub("-", ".", names) # Names and conditions do not tolerate "-"
conditions <- unlist(tstrsplit(args[5], ","))
conditions <- gsub("-", ".", conditions) # Names and conditions do not tolerate "-"
controls <- unlist(tstrsplit(args[6], ","))
controls <- gsub("-", ".", controls) # Names and conditions do not tolerate "-"
dds_output_folder <- args[7]
FC_output_folder <- args[8]
PDF_output_folder <- args[9]
experiment <- args[10]
feature <- args[11]
norm <- args[12]

# Import data ----
dat <- lapply(counts, fread)
names(dat) <- names
dat <- rbindlist(dat, idcol = "condition")
dat[, condition:= factor(condition, unique(condition))]
DF <- dcast(dat, ID~condition, value.var = "count")
DF <- data.frame(DF[, -1], row.names = DF$ID)

# Remove low count reads ----
DF <- DF[rowSums(DF >= 3) >= 2,]

# SampleTable ----
sampleTable <- data.frame(condition = conditions,
                          row.names = names)

# Ref genome read counts ----
ref <- lapply(refStats, fread)
names(ref) <- names
ref <- rbindlist(ref, idcol = "sample")

# Spike-in read counts ----
spike <- lapply(spikeStats, fread)
names(spike) <- names
spike <- rbindlist(spike, idcol = "sample")

# Merge statistics ----
stats <- merge(ref, spike, by= "sample", suffixes= c("", "_spikeIn"), sort= FALSE)
stats[, cdition:= conditions]

# DESeq2 analysis ----
print(paste("Start", norm, "normalization"))
dds <- DESeqDataSetFromMatrix(countData = DF,
                              colData = sampleTable,
                              design = ~ condition)

# SizeFactors ----
if(norm=="libSize") {
  sizeFactors(dds) <- stats[, umi_counts/median(umi_counts)]
} else if(norm=="spikeIn") {
  sizeFactors(dds) <- stats[, umi_counts_spikeIn/median(umi_counts_spikeIn)]
} else if(norm!="default") {
  stop("normalization should be one of 'default', 'libSize' or 'spikeIn'")
}

# Compute model and save object ----
dds <- DESeq(dds)
saveRDS(dds,
        paste0(dds_output_folder, "/", experiment, "_", feature, "_", norm, "_norm_DESeq2.dds"))

# Open pdf to save MA plot ----
outputPdf <- paste0(PDF_output_folder, "/", experiment, "_", feature, "_", norm, "_norm_MAplots.pdf")
pdf(outputPdf, 4, 3)
par(mai= c(.9,1.5,.9,1.3),
    cex.axis= 6/12,
    cex.lab= 7/12,
    las= 1,
    tcl= -0.1,
    mgp= c(.8, 0.25, 0),
    font.main= 1,
    cex.main= 9/12)

# Initiate list to save FC ----
FC <- list()

# For each control and each condition, compute FC and print MA plot ----
for(ctl in unique(controls))
{
  for(cdition in setdiff(unique(conditions), ctl))
  {
    # Compute FC table
    res <- results(dds,
                   contrast = c("condition", cdition, ctl))
    res <- as.data.frame(res)
    res <- as.data.table(res, keep.rownames = "gene_id")
    # Compute diff column
    res[, diff:= fcase(padj<0.05 & log2FoldChange>log2(1.5), "Up-regulated",
                       padj<0.05 & log2FoldChange<(-log2(1.5)), "Down-regulated",
                       default = "Unaffected")]
    # Add condition and control columns
    res[, condition:= cdition]
    res[, control:= ctl]
    # Add to FC list
    FC[[paste0(cdition, "_vs_", ctl)]] <- data.table::copy(res)

    # MA plot colors
    res[, col:= {
      fcase(diff=="Up-regulated", "tomato",
            diff=="Down-regulated", "cornflowerblue",
            default = "lightgrey")
    }]

    # Plot unaffected genes first
    res <- res[order(col=="lightgrey", decreasing = TRUE)]

    # Compute limits and pch
    lims <- quantile(res$log2FoldChange, c(0.001, 0.999))
    res[, pch:= ifelse(between(log2FoldChange, lims[1], lims[2]), 16, 17)]

    # Clip outliers
    res[log2FoldChange<lims[1], log2FoldChange:= lims[1]]
    res[log2FoldChange>lims[2], log2FoldChange:= lims[2]]

    # Plot
    res[, {
      plot(x = log10(baseMean),
           y = log2FoldChange,
           col= adjustcolor(col, .5),
           pch= pch,
           ylab= "Fold change (log2)",
           frame= F,
           xaxt= "n",
           cex= .5,
           main= paste(cdition, "vs.", ctl))
      axis(1, padj= -1.45)
      abline(h= 0, lty= 3)

      # Legend
      nUp <- sum(diff=="Up-regulated")
      nUp <- formatC(nUp, big.mark = ",")
      nDown <- sum(diff=="Down-regulated")
      nDown <- formatC(nDown, big.mark = ",")
      legend(par("usr")[2],
             par("usr")[4],
             col= adjustcolor(c("tomato", "cornflowerblue"), 0.5),
             legend= c(paste0("Up-regulated (", nUp, ")"),
                       paste0("Down-regulated (", nDown, ")")),
             pch= 16,
             bty= "n",
             cex= 6/12,
             border= NA,
             xpd= NA)
    }]
  }
}
dev.off()

# Save FC file ----
FC <- rbindlist(FC)
setcolorder(FC, c("condition", "control"))
outputFile <- paste0(FC_output_folder, "/", experiment, "_", feature, "_", norm, "_norm_DESeq2_FC.txt")
fwrite(FC,
       outputFile,
       sep="\t",
       na = NA)
