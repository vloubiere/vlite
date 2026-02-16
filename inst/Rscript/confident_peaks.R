#!/usr/bin/env Rscript

# Check whether required arguments are provided ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage:
       [required] 1/ replicates .narrowPeak files (comma separated) \n
       [required] 2/ merged .narrowPeak (or .broadPeak) file \n
       [required] 3/ output_file.narrowPeak (or .broadPeak) \n
       [required] 4/ output_file.pdf path \n")
}

# Load packges ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(devtools, warn.conflicts = FALSE))
devtools::load_all("/zssd/scratch/vincent.loubiere/vlite/")

# Parse arguments ----
rep.files <- unlist(tstrsplit(args[1], ","))
merged.file <- args[2]
output.file <- args[3]
pdf.file <- args[4]

# Import replicates peaks ----
reps <- lapply(rep.files, function(x) {
  .c <- fread(x,
              select = c(1,2,3,9),
              col.names = c("seqnames", "start", "end", "qValue"))
  .c[, start:= start+1] # 0 base to 1 base
  return(.c)
})

# Import Merged peaks ----
merge <- fread(merged.file)
col.names <- c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
setnames(merge,
         col.names[1:ncol(merge)])
merge[, start:= start+1] # 0 base to 1 base

# Initiate plot ----
pdf(pdf.file, 8, 4)
vl_par(mfrow= c(1,2))

# Check ----
if(nrow(merge)) {
  # Compute mean qValue per replicate ----
  for(i in seq(reps)) {
    # Name
    .n <- paste0("rep.", i)
    # Data
    .c <- reps[[i]][merge, mean(qValue), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1
    merge[, (.n):= .c]
    if(i>1) {
      x <- merge$rep.1
      x[is.na(x)] <- 0
      y <- merge[[.n]]
      y[is.na(y)] <- 0
      
      vl_plot(
        x,
        y,
        col= adjustcolor(ifelse(x==0 | y==0, "grey", "red"), .3),
        xlab= "-log10(qValue) rep 1",
        ylab= paste0("-log10(qValue) rep ", i),
        cex= .5
      )
      addPcc(cor(x,y))
      title(main= gsub(".pdf$", "", basename(pdf.file)), cex.main= 1, cex.font= 1, line = 2)
    }
  }
  
  # Upset plot ----
  pl <- melt(merge,
             id.vars = "name",
             measure.vars = patterns("^rep."))[!is.na(value)]
  par(mai= c(1.5, 1.5, .9, .5))
  upsetPlot(split(pl$name, pl$variable))
} else {
  plot.new()
  text(.5, .5, "No peaks")
}
dev.off()

# Select confident peaks ----
merge[, rm:= apply(.SD, 1, anyNA), .SDcols= patterns("^rep.")]
final <- merge[(!rm), !c("rm", grep("^rep.", names(merge), value= T)), with= FALSE]

# Save ----
final[, start:= start-1] # 0-based bed
fwrite(final,
       output.file,
       col.names = F,
       quote= F,
       sep= "\t",
       na = ".")
