#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript script.R {replicates .narrowPeak files (comma separated)} {merged .narrowPeak (or .broadPeak) file} {output_file.narrowPeak (or .broadPeak)}")
}

# Load packges ----
suppressMessages(library(data.table, warn.conflicts = FALSE))

# Parse arguments ----
rep.files <- strsplit(args[[1]], ",")[[1]]
merged.file <- args[[2]]
output.file <- args[[3]]

# Import replicates peaks ----
reps <- lapply(rep.files, function(x) {
  .c <- fread(x,
              select = 1:3,
              col.names = c("seqnames", "start", "end"))
  .c[, start:= start+1] # 0 base to 1 base
  return(.c)
})

# Import Merged peaks ----
merge <- fread(merged.file)
col.names <- c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
setnames(merge,
         col.names[1:ncol(merge)])
merge[, start:= start+1] # 0 base to 1 base
conf <- data.table::copy(merge)

# Extract confident peaks found in all replicates and save ----
for(i in seq(reps)) {
  sel <- reps[[i]][conf, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N>0
  conf <- conf[sel]
}

# Save ----
conf[, start:= start-1] # 0-based bed
fwrite(conf,
       output.file,
       col.names = F,
       quote= F,
       sep= "\t",
       na = ".")
