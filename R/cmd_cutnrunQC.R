#' Generate commands for peaks QC
#'
#' @description
#' Creates commands to do QC of individual CUTNRUN/ChIP-seq/ATAC-Seq peaks.
#'
#' @param peaks.file Path to the sample's .narrowpeak or .broadPeak file.
#' @param bw.sample Path to sample bigwig.
#' @param bw.input Optional path to the input bigwig file that was used for peak calling.
#' @param genome A BSgenome name (e.g. "mm10", "dm6"...).
#' @param output.prefix Prefix for the output files. If not specified, will be generated from the bw.sample file name.
#' @param pdf.output.folder Output folder where the pdf file should be saved. Default= "db/QC/CUTNRUN/".
#' @param Rpath Path to the Rscript binary. Default= "Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file labels (e.g: "peaks", "bedgraph").
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run the peak calling pipeline.
#' - `job.name`: Default name for the job = "MACS2".
#'
#' @examples
#'
#' @export
cmd_cutnrunQC <- function(peaks.file,
                          bw.sample,
                          bw.input= NULL,
                          genome,
                          output.prefix,
                          pdf.output.folder= "db/QC/CUTNRUN/",
                          Rpath = "Rscript")
{
  # Check (!Do not check if bam file(s) exist to allow wrapping!) ----
  if(length(peaks.file)>1)
    stop("peaks.file file should be unique.")
  if(!grepl(".narrowPeak$|.broadPeak$", peaks.file))
    stop("Peaks file should be in .narrowPeak or .broadPeak format.")
  if(length(bw.sample)>1)
    stop("bw.sample file should be unique.")
  if(length(bw.input)>1)
    stop("bw.input file should be unique.")
  if(missing(output.prefix))
    output.prefix <- gsub(".bw$", "", basename(bw.sample))

  # Output files paths ----
  pdf.file <- file.path(pdf.output.folder, paste0(output.prefix, "_QC.pdf"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "cutnrunQC.R", package = "vlite"),
    peaks.file, # Peaks file
    pdf.file, # Output pdf file
    genome, # BSgenome name
    bw.sample, # Sample bigwig file
    bw.input # Input bigwig file
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "peaks.QC.pdf",
                    path= pdf.file,
                    cmd= cmd)

  # Return ----
  cmd$job.name <- "peaksQC"
  return(cmd)
}
