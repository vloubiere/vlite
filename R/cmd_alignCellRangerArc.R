#' Generate Cell Ranger ARC Alignment Command for Multiome Data
#'
#' @description
#' Constructs a shell command to run 10x Genomics Cell Ranger ARC for paired
#' multiome (GEX and ATAC) quantification. This function is intended for
#' real data and does not create dummy files.
#'
#' @param gex.dir Path to the directory containing GEX FASTQ files.
#' @param gex.sample.name The sample name prefix for GEX FASTQs.
#' cellranger-arc will only search for filenames matching sampleName_S*_L.*.
#' @param atac.dir Path to the directory containing ATAC FASTQ files.
#' @param atac.sample.name The sample name prefix for ATAC FASTQs.
#' cellranger-arc will only search for filenames matching sampleName_S*_L.*.
#' @param index Path to the Cell Ranger ARC reference.
#' @param output.prefix Prefix for output files and the run ID.
#' @param output.folder Directory where output files will be saved. Default= 'db/multiome_10X/'.
#' @param cores Number of CPU cores to use. Default= 10.
#' @param mem Memory in GB to allocate. Default= 64.
#'
#' @export
cmd_alignCellRangerArc <- function(gex.dir,
                                   gex.sample.name,
                                   atac.dir,
                                   atac.sample.name,
                                   index,
                                   output.prefix,
                                   output.folder = "db/multiome_10X/",
                                   cores = 10,
                                   mem = 64)
{
  # Checks ----
  gex.dir <- unique(gex.dir)
  gex.sample.name <- unique(gex.sample.name)
  atac.dir <- unique(atac.dir)
  atac.sample.name <- unique(atac.sample.name)
  if(length(gex.dir)>1 | length(gex.sample.name)>1 | length(atac.dir)>1 | length(atac.sample.name)>1)
    stop("gex.dir, gex.sample.name, atac.dir or atac.sample.name should all be unique.")
  # A GEX library should have at least R1 and R2 reads.
  if(length(list.files(gex.dir, "fastq.gz$")) < 2)
    stop("Less than 2 fastq.gz files in the provided gex.dir.")
  # A Multiome ATAC library has R1, R2, and I1 reads.
  if(length(list.files(atac.dir, "fastq.gz$")) < 3)
    stop("Less than 3 fastq.gz files in the provided atac.dir. A multiome ATAC run requires R1, R2, and I1 files.")

  # Create output folder to save the library files ----
  if (!dir.exists(output.folder)) {
    dir.create(output.folder, recursive = TRUE)
  }

  # Create the libraries CSV with real GEX and real ATAC data ----
  libraries_df <- data.frame(
    fastqs = c(normalizePath(gex.dir), normalizePath(atac.dir)), # FIXED: Was using dummy.atac.dir
    sample = c(gex.sample.name, atac.sample.name),
    library_type = c("Gene Expression", "Chromatin Accessibility")
  )
  # Save
  libraries.csv.path <- file.path(output.folder, paste0(output.prefix, "_libraries.csv"))
  write.csv(libraries_df,
            file = libraries.csv.path,
            row.names = FALSE,
            quote = FALSE)

  # Define paths for both GEX and ATAC output BAM files ----
  gex.bam.path <- file.path(output.folder, output.prefix, "outs", "gex_possorted_bam.bam")
  atac.bam.path <- file.path(output.folder, output.prefix, "outs", "atac_possorted_bam.bam")

  # Construct the cellranger-arc command ----
  cmd <- paste0(
    "module load build-env/f2022; ",
    "module load cellranger-arc/2.0.2; ",
    "cd ", normalizePath(output.folder), "; ",
    "cellranger-arc count --id=", output.prefix,
    " --reference=", normalizePath(index),
    " --libraries=", normalizePath(libraries.csv.path),
    " --localcores=", cores,
    " --localmem=", mem
  )

  # Wrap commands output, creating one row for each BAM file ----
  cmd <- data.table::data.table(
    file.type = c("lib.csv", "bam.gex", "bam.atac"),
    path = c(libraries.csv.path, gex.bam.path, atac.bam.path),
    cmd = cmd,
    cores = cores,
    mem = mem,
    job.name = "cellRangerArc"
  )

  # Return ----
  return(cmd)
}
