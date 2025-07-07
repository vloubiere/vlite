#' Convert BED to BigWig Format
#'
#' @description
#' Creates shell commands to convert a BED file to a CPM-normalized bigwig file.
#'
#' @param bed Path to the input BED file. Only a single BED file is allowed.
#' @param genome A BS genome name use to retrieve seqnames.
#' @param bed.subset Optional path to a bed file. If provided, only reads overlapping these regions on the same
#' strand are used. Default= NULL (no filtering).
#' @param output.prefix Prefix for the output BigWig file. If not provided, it is derived from the input bed filename.
#' @param bw.output.folder Directory for the BigWig file. Default= "db/bw/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file label ("bw").
#' - `path`: Path to the BigWig file.
#' - `cmd`: Shell command to run the BED to BigWig conversion.
#' - `job.name`: Default name for the job = "bedToBw".
#'
#' @examples
#' # Convert a BED file to BigWig format
#' cmd <- cmd_bedToBigwig(
#'   bed = "/data/bed/sample.bed",
#'   genome = "mm11"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_bedToBigwig <- function(bed,
                            genome,
                            bed.subset= NULL,
                            output.prefix= NULL,
                            bw.output.folder= "db/bw/",
                            Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if bed file exists to allow wrapping!) ----
  if(length(bed)!=1)
    stop("A unique bed file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bed$", "", basename(bed))

  # Output file ----
  bw.file <- file.path(bw.output.folder, paste0(output.prefix, "_", genome, ".bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "bed_to_bigwig.R", package = "vlite"),
    bed,
    genome,
    bw.file,
    bed.subset
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= bw.file,
                    cmd= cmd,
                    job.name= "bedToBw")

  # Return ----
  return(cmd)
}
