#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Please specify:\n
       [required] 1/ A comma-separated list of fq1 files \n
       [required] 2/ A comma-separated list of fq2 files \n
       [required] 3/ subreadr index prefix \n
       [required] 4/ Output bam path \n")
}

# Load packges ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(Rsubread, warn.conflicts = FALSE))

# Parse arguments ----
fq1 <- unlist(tstrsplit(args[1], ","))
fq2 <- unlist(tstrsplit(args[2], ","))
if(fq2=="")
  fq2 <- NULL
idx <- args[3]
bam <- args[4]

# Catenate fq files if necessary
if(length(fq1)>1)
{
  temp_file <- tempfile(fileext = ".txt")
  system(paste("cat", paste(fq1, collapse = " "), ">", temp_file))
  fq1 <- temp_file
}
if(!is.null(fq2) && length(fq2)>1)
{
  temp_file <- tempfile(fileext = ".txt")
  system(paste("cat", paste(fq2, collapse = " "), ">", temp_file))
  fq2 <- temp_file
}

# Align ----
Rsubread::align(index= idx,
                readfile1= fq1,
                readfile2= fq2,
                type= "rna",
                input_format = "gzFASTQ",
                output_format = "BAM",
                maxMismatches = 3,
                nthreads = data.table::getDTthreads()-1,
                unique = T,
                output_file= bam)
