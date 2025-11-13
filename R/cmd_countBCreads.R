#' Generate Commands for Counting Barcode Reads
#'
#' @description
#' Creates shell commands to count barcode reads from a BAM file using a predefined barcode dictionary.
#' Outputs a counts file.
#'
#' @param bam Path to the input BAM file. Only a single BAM file is allowed.
#' @param output.prefix Prefix for the output file. If not provided, it is derived from the input BAM filename.
#' @param counts.output.folder Directory for the counts file. Default= "db/counts/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file label ("counts.BC").
#' - `path`: Path to the output counts file.
#' - `cmd`: Shell command to run the barcode counting pipeline.
#' - `job.name`: Default name for the job = "countBCreads".
#'
#' @examples
#' # Count barcodes in a BAM file
#' cmd <- cmd_countBCreads(
#'   bam = "/data/bam/sample.bam"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_countBCreads <- function(bam,
                             output.prefix= NULL,
                             counts.output.folder= "db/counts/",
                             Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if bam file exists to allow wrapping!) ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))

  # Output file ----
  counts.file <- file.path(counts.output.folder, paste0(output.prefix, "_counts.txt"))

  # Command ----
  cmd <- paste(Rpath,
               system.file("Rscript", "BC_counts.R", package = "vlite"),
               bam,
               counts.file)

  # Wrap commands output ----
  cmd <- data.table(file.type= "counts.BC",
                    path= counts.file,
                    cmd= cmd,
                    job.name= "countBCreads")

  # Return ----
  return(cmd)
}
