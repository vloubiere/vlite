#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% 6:7) {
  stop("Please specify:\n
       [required] 1/ Input bam path \n
       [required] 2/ Layout ('PAIRED' or 'SINGLE')
       [required] 3/ Path to a gtf or a saf annotation file (ideally the one used to generate the subreadr index) \n
       [required] 4/ Output statistics file \n
       [required] 5/ Output count file \n
       [required] 6/ If a read overlaps more than one feature, should be assigned to all overlapping features? \n
       [required] 7/ [Optional] The name of the an extra column to be reported in the counts files (e.g. 'gene_symbol') \n")
}

# Load packges ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(Rsubread, warn.conflicts = FALSE))

# Parse arguments ----
bam <- args[1]
layout <- args[2]
if(!layout %in% c("PAIRED", "SINGLE"))
  stop("layout shoud be one of 'PAIRED' or 'SINGLE'")
annot.ext <- args[3]
isGTFAnnotationFile <- grepl(".gtf$|.gtf.gz$", annot.ext)
# Check format
if(!isGTFAnnotationFile) {
  if(!grepl(".saf$", annot.ext))
    stop("annot.ext should be a '.gtf', '.gtf.gz' or a '.saf' file")
  check <- names(fread(annot.ext, nrows = 0))
  if(!all(c("GeneID", "Chr", "Start", "End", "Strand") %in% check))
    stop(".saf file should contain GeneID, Chr, Start, End and Strand columns.")
}
stats_file <- args[4]
count_file <- args[5]
allowMultiOverlap <- as.logical(args[6])
GTF.attrType.extra <- if(length(args)==7)
  args[7] else
    NULL

# Count ----
.c <- Rsubread::featureCounts(files= bam,
                              annot.ext= annot.ext,
                              isGTFAnnotationFile = isGTFAnnotationFile,
                              GTF.featureType= "exon",
                              GTF.attrType= "gene_id",
                              GTF.attrType.extra= if(isGTFAnnotationFile) GTF.attrType.extra else NULL,
                              allowMultiOverlap = allowMultiOverlap,
                              isPairedEnd = layout=="PAIRED",
                              nthreads = max(c(1, data.table::getDTthreads()-1)))

# Save statistics ----
count_stats <- .c$stat
setnames(count_stats, c("feature", "count"))
fwrite(count_stats,
       stats_file,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote= F)

# Assemble full gene id ----
gene_id <- as.data.table(.c$annotation)
if(isGTFAnnotationFile && !is.null(GTF.attrType.extra)) {
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
