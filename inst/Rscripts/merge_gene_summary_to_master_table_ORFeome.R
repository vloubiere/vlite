#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Please specify:\n
       [required] 1/ Gene summary output file from mageck \n
       [required] 2/ Path to master_table \n
       [required] 3/ sort. can be one of 'pos' or 'neg'. \n
       [required] 4/ Output file path \n")
}

require(data.table)

# Example ----
# gene_summary <- "db/FC_tables/ORFeome/apoptosisFasL_A549/apoptosisFasL_A549_FasL.pos.gene_summary.txt"
# output_file <- "db/FC_tables/ORFeome/apoptosisFasL_A549/apoptosisFasL_A549_FasL.pos.gene_summary_master_posFC.txt"
# sort <- "pos"
# master_table <- "/groups/stark/pachano/projects/eORFeome/Rdata/Master_eORFeome_Mar24.csv"

# parse ----
gene_summary <- args[1]
master_table <- args[2]
sort <- args[3]
output_file <- args[4]

# Checks ----
if(!grepl(".txt$", gene_summary))
  stop("gene_summary should be a path to a .txt file outputed by MAGeCK")
if(!sort %in% c("pos", "neg"))
  stop("sort should be one of pos or neg!")
if(!grepl(".csv$", master_table))
  stop("master should be a path to a .csv file")

# Import master ----
master <- fread(master_table)
master <- as.data.table(master)
master[, id:= as.character(id)]
# Columns of interest
# master.columns <- c("id", "Effector class", "RefSeq protein", "UniProt", "Gene",
#                     "Description", "Species", "DNA sequence", "size", "Amino acid sequence",
#                     "Size (aa)", "Uniprot", "Protein names", "Function")
master.columns <- c("id", "Class", "sp_abb", "Species", "seq",
                    "size", "Uniprot", "Family", "Group", "Group",
                    "genome")
master <- master[, (master.columns), with= F]

# Import FC table ----
FC <- fread(gene_summary)
if(!"id" %in% names(FC))
  stop("gene_summary should contain an id column")
FC[, first_id:= tstrsplit(id, "_", keep= 1)]

# Merge tables ----
# Cols of interest
cols <- c("first_id", "id", "num")
cols <- if(sort=="pos")
  c(cols, "pos|rank", "pos|goodsgrna", "pos|score", "pos|lfc", "pos|fdr") else
    c(cols, "neg|rank", "neg|goodsgrna", "neg|score", "neg|lfc", "neg|fdr")
# merge
res <- merge(FC[, c("first_id", "id", "num", "pos|rank", "pos|goodsgrna", "pos|score", "pos|lfc", "pos|fdr")],
             master,
             by.x= "first_id",
             by.y= "id",
             all.x= TRUE)
# simplify names
if(sort=="pos")
  setnames(res, function(x) gsub("pos|", "", fixed = TRUE, x)) else
    setnames(res, function(x) gsub("neg|", "", fixed = TRUE, x))

# Save ----
fwrite(res,
       output_file,
       col.names = T,
       row.names = F,
       na= NA,
       sep = "\t",
       quote=F)
