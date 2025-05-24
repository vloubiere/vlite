#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Please specify:\n
       [required] 1/ Merged FC table (.txt). \n
       [required] 2/ pdf output file (.pdf). \n")
}

require(ggplot2)
require(ggrepel)
require(data.table)

# Import data ----
FC_table <- args[1]
pdf <- args[2]

# Checks ----
stopifnot(grepl(".txt$", FC_table))
stopifnot(grepl(".pdf$", pdf))

# Import FC data
FC <- fread(FC_table)

# Prepare for plotting ----
FC[, logFDR:= -log10(fdr)]
yMax <- quantile(FC$lfc, .999, na.rm= T)
FC[, shape:= ifelse(logFDR > yMax, "triangle", "circle")]
FC[logFDR>yMax, logFDR:= yMax]
FC[, col:= ifelse(hit, "Hit", "None")]
plMageck <- ggplot(FC, aes(x = lfc, y = logFDR)) +
  geom_point(aes(color = col, shape = shape)) +
  geom_text_repel(data = FC[(col!="None")],
                  aes(label = id, col= col),
                  max.overlaps = Inf,
                  size = 2) +
  theme_minimal() +
  labs(title = gsub("_FC_MAGeCK.txt$", "", basename(FC_table)),
       x = "Fold Change (log)",
       y = "FDR (-log10)") +
  scale_color_manual(values = c("Hit" = "red", "None" = "lightgrey"), name = "Significant")+
  ylim(0, yMax) +
  guides(shape = "none")  # Remove shape legend

# Print pdf
pdf(pdf)
plot(plMageck)
dev.off()
