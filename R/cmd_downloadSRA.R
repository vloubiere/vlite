#' Download and Process SRA Files
#'
#' @param SRR SRA accession number
#' @param layout Either "SINGLE" or "PAIRED" for read layout
#' @param output.prefix Output file prefix (default= SRR number).
#' @param fq.output.folder Directory for output FASTQ files.
#' @param gzip Logical, whether to gzip the output files.
#' @param tmp.folder Directory for temporary files.
#' @param threads Number of threads (default: 4)
#' @param mem Memory in MB (default: 2048)
#'
#' @return data.table with columns:
#'   \item{file.type}{Type of output file (fq1 or fq2)}
#'   \item{path}{Full path to output file}
#'   \item{cmd}{Command to execute}
#'
#' @examples
#' cmd <- cmd_downloadSRA("SRR123456",
#'                        layout="PAIRED",
#'                        fq.output.folder="./fastq",
#'                        tmp.folder="./tmp")
#' # Submit command
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_downloadSRA <- function(SRR,
                            layout,
                            output.prefix= SRR,
                            fq.output.folder= "db/fq/",
                            gzip= TRUE,
                            tmp.folder= paste0(fq.output.folder, "/tmp/"),
                            threads= 4,
                            mem= 2048)
{
  # Checks ----
  if(length(SRR)!=1)
    stop("A unique SRR id should be provided.")
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")

  # Downloaded files ----
  if(layout=="PAIRED") {
    fq1 <- file.path(fq.output.folder, paste0(SRR, "_1.fastq"))
    fq2 <- file.path(fq.output.folder, paste0(SRR, "_2.fastq"))
  } else if(layout=="SINGLE") {
    fq1 <- file.path(fq.output.folder, paste0(SRR, ".fastq"))
  }

  # Renamed files ----
  if(layout=="PAIRED") {
    fq1.new <- file.path(fq.output.folder, paste0(output.prefix, "_1.fq"))
    fq2.new <- file.path(fq.output.folder, paste0(output.prefix, "_2.fq"))
  } else if(layout=="SINGLE") {
    fq1.new <- file.path(fq.output.folder, paste0(output.prefix, ".fq"))
  }

  # Download command ----
  cmd <- "module load build-env/f2022; module load sra-toolkit/3.1.1-centos_linux64; fasterq-dump"
  if(layout=="PAIRED")
    cmd <- paste(cmd, "--split-files")
  cmd <- paste(cmd, "--threads", threads, "--mem", mem, "--temp", tmp.folder, "--outdir", fq.output.folder, SRR)

  # Automatically gzip the output files ----
  cmd <- paste(cmd, "; mv", fq1, fq1.new)
  if(layout=="PAIRED")
    cmd <- paste(cmd, "; mv", fq2, fq2.new)

  # Automatically gzip the output files ----
  if(gzip) {
    cmd <- paste(cmd, "; gzip", file.path(fq.output.folder, paste0(output.prefix, "*.fq")))
    fq1.new <- paste0(fq1.new, ".gz")
    if(layout=="PAIRED")
      fq2.new <- paste0(fq2.new, ".gz")
  }

  # Wrap commands output ----
  cmd <- data.table(file.type= if(layout=="PAIRED") c("fq1", "fq2") else "fq1",
                    path= if(layout=="PAIRED") c(fq1.new, fq2.new) else fq1.new,
                    cmd= cmd)

  # Return
  return(cmd)
}
