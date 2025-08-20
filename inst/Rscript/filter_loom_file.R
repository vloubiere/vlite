#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Test if there are 5 args: if not, return an error ----
if (length(args) != 4) {
  stop(
    "Please specify:\n
       [required] 1/ The path to a Seurat object \n
       [required] 2/ The path to a loom file \n
       [required] 3/ Condition name \n
       [required] 4/ Output folder \n"
  )
}

# Load libraries ----
suppressMessages(library(Seurat, warn.conflicts = FALSE))
suppressMessages(library(SeuratWrappers, warn.conflicts = FALSE))
suppressMessages(library(velocyto.R, warn.conflicts = FALSE))
suppressMessages(library(loomR, warn.conflicts = FALSE))
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Examples
# args <- c("db/seurat_individual_samples/final/PH18_final.rds",
#           "db/snRNASeq_10X/PH18/PH18.loom",
#           "PH18",
#           "db/velocity/")

# Parse arguments
seurat.object <- args[1]
loom.file <- args[2]
cdition <- args[3]
output.folder <- file.path(args[4], cdition)
dir.create(output.folder, showWarnings = F, recursive = T)

# Output files
output.file <- file.path(output.folder, paste0(cdition, "_filtered.loom"))
if(file.exists(output.file))
  file.remove(output.file)

# Load Seurat object
seurat_v5 <- readRDS(seurat.object)

# Load loom file
ldat <- read.loom.matrices(loom.file)

# Subset matrices cells and genes present in seurat object
filtered <- lapply(ldat, function(x) {
  colnames(x) <- gsub(paste0(cdition, ":(.*)x"), "\\1-1", colnames(x))
  return(x[rownames(seurat_v5), colnames(seurat_v5)])
})

# Create loom with spliced as primary 'matrix'
lf <- loomR::create(
  filename = output.file,
  data = filtered$spliced,
  cell.attrs = list(cell_cluster= Idents(seurat_v5),
                    UMAP= as.matrix(seurat_v5@reductions$umap@cell.embeddings)),
  layers = list(
    spliced= filtered$spliced,
    unspliced= filtered$unspliced
  ),
  overwrite = TRUE,
  verbose = TRUE
)
lf$close_all()

# Save pca as csv
fwrite(
  as.data.table(seurat_v5@reductions$pca@cell.embeddings, keep.rownames = "cell"),
  file.path(output.folder, "pca_from_R.csv"),  # columns: cell, PC1..PCk  
)

# Save umap as csv
fwrite(
  as.data.table(seurat_v5@reductions$umap@cell.embeddings, keep.rownames = "cell"),
  file.path(output.folder, "umap_from_R.csv"),  # columns: cell, PC1..PCk  
)
