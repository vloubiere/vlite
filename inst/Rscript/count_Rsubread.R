#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% c(5,6)) {
  stop("Please specify:\n
       [required] 1/ Input bam path \n
       [required] 2/ Layout ('PAIRED' or 'SINGLE')
       [required] 3/ Path to a gtf file (ideally the one used to generate the subreadr index) \n
       [required] 4/ Output statistics file \n
       [required] 5/ Output count file \n
       [required] 6/ [Optional] The name of the an extra column to be reported in the counts files (e.g. 'gene_symbol') \n")
}

# Load packges ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(Rsubread, warn.conflicts = FALSE))

# Parse arguments ----
bam <- args[1]
layout <- args[2]
if(!layout %in% c("PAIRED", "SINGLE"))
  stop("layout shoud be one of 'PAIRED' or 'SINGLE'")
gtf <- args[3]
stats_file <- args[4]
count_file <- args[5]
GTF.attrType.extra <- if(length(args)==6)
  args[6] else
    NULL

# Count ----
.c <- Rsubread::featureCounts(bam,
                              annot.ext= gtf,
                              isGTFAnnotationFile = TRUE,
                              GTF.featureType= "exon",
                              GTF.attrType= "gene_id",
                              GTF.attrType.extra= GTF.attrType.extra,
                              isPairedEnd = layout=="PAIRED",
                              nthreads = data.table::getDTthreads()-1)

# Save statistics ----
count_stats <- .c$stat
setnames(count_stats, c("feature", "count"))
fwrite(count_stats,
       stats_file,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote= F)

# Assmeble full gene id ----
gene_id <- as.data.table(.c$annotation)
if(!is.null(GTF.attrType.extra)) {
  gene_id[, GeneID:= paste0(GeneID, "__", gene_id[[GTF.attrType.extra]])]
}
gene_id[, GeneID:= paste0(GeneID, "__", Length)]

# Save count table ----
counts <- data.table(gene_id= gene_id$GeneID,
                     count= .c[[1]][,1])
fwrite(counts,
       count_file,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote= F)
