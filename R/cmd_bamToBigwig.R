#' Convert BAM to BigWig Format
#'
#' @description
#' Creates shell commands to convert a BAM file to BigWig format. Supports single-end and paired-end sequencing data,
#' with optional fragment extension for paired-end reads.
#'
#' @param bam Path to the input BAM file. Only a single BAM file is allowed.
#' @param layout Sequencing layout, either "SINGLE" or "PAIRED".
#' @param output.prefix Prefix for the output BigWig file. If not provided, it is derived from the input BAM filename.
#' @param extend.PE.fragments Logical. Whether to extend paired-end fragments. Default= FALSE.
#' @param extsize Numeric. Read extension size. Default= 0 (no extension).
#' @param bw.output.folder Directory for the BigWig file. Default= "db/bw/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of CPU cores to use. Default= 6.
#'
#' @return A data.table with:
#' - `file.type`: Output file label ("bw").
#' - `path`: Path to the BigWig file.
#' - `cmd`: Shell command to run the BAM to BigWig conversion.
#' - `job.name`: Default name for the job = "bamToBw".
#'
#' @examples
#' # Convert a BAM file to BigWig format for single-end data
#' cmd <- cmd_bamToBigwig(
#'   bam = "/data/bam/sample.bam",
#'   layout = "SINGLE"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Convert a BAM file to BigWig format for paired-end data with fragment extension
#' cmd <- cmd_bamToBigwig(
#'   bam = "/data/bam/sample.bam",
#'   layout = "PAIRED",
#'   extend.PE.fragments = TRUE
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_bamToBigwig <- function(bam,
                            layout,
                            output.prefix= NULL,
                            extend.PE.fragments= FALSE,
                            extsize= 0,
                            bw.output.folder= "db/bw/",
                            Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                            cores= cores)
{
  # Check (!Do not check if bam file exists to allow wrapping!) ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(extsize<0)
    stop("extsize < 0 not supported")

  # Output file ----
  bw.file <- file.path(bw.output.folder, paste0(output.prefix, ".bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "bam_to_bigwig.R", package = "vlite"),
    bam,
    layout,
    extend.PE.fragments,
    bw.file,
    extsize
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= bw.file,
                    cmd= cmd,
                    cores= cores,
                    job.name= "bamToBw")

  # Return ----
  return(cmd)
}
