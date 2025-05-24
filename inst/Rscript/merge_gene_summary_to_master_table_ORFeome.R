#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 6) {
  stop("Please specify:\n
       [required] 1/ Gene summary output file from MAGeCK (.txt). \n
       [required] 2/ sort. can be one of 'pos' or 'neg'. \n
       [required] 3/ logFC cutoff to define hits. \n
       [required] 4/ FDR cutoff to define hits. \n
       [required] 5/ Path to .rds master_table. Should contain at least and 'id' column (ORF ids). \n
       [required] 6/ Output file path (.txt). \n")
}

require(data.table)

# Example ----
# gene_summary <- "db/FC_tables/ORFeome/cutoffs_s3_i0_t3_p1_param/NFkB_Jurkat.pos.gene_summary.txt"
# output_file <- "db/FC_tables/ORFeome/apoptosisFasL_A549/apoptosisFasL_A549_FasL.pos.gene_summary_master_posFC.txt"
# sort <- "pos"
# master_table <- "/groups/stark/pachano/projects/eORFeome/Rdata/Master_eORFeome_Mar24.csv"

# parse ----
gene_summary <- args[1]
sort <- args[2]
logFCcutoff <- as.numeric(args[3])
FDRcutoff <- as.numeric(args[4])
master_table <- args[5]
output_file <- args[6]

# Checks ----
if(!grepl(".txt$", gene_summary))
  stop("gene_summary should be a path to a .txt file outputed by MAGeCK")
if(!sort %in% c("pos", "neg"))
  stop("sort should be one of pos or neg!")
if(!grepl(".rds$", master_table))
  stop("master should be a path to a .rds file")

# Import FC table ----
FC <- fread(gene_summary)

# Filter columns ----
cols <- if(sort=="pos") {
  c("id", "pos|lfc", "pos|fdr", "pos|rank", "pos|score", "num", "pos|goodsgrna")
} else if(sort=="neg") {
  c("id", "neg|lfc", "neg|fdr", "neg|rank", "neg|score", "num", "neg|goodsgrna")
}
FC <- FC[, cols, with= FALSE]
setnames(FC, c("id", "lfc", "fdr", "rank", "score", "num", "goodsgrna"))
if(sort=="neg")
  FC[, lfc:= -lfc]

# Define hits ----
FC[, hit:= lfc >= logFCcutoff & fdr <= FDRcutoff]

# Import master ----
master <- readRDS(master_table)

# Merge tables and save ----
FC[, id:= as.character(id)]
res <- merge(master,
             FC,
             by= "id",
             all.x= TRUE)

# Save ----
fwrite(res,
       output_file,
       col.names = T,
       row.names = F,
       na= NA,
       sep = "\t",
       quote=F)
