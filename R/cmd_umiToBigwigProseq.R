#' Convert UMI Counts to BigWig Format for PRO-Seq Data
#'
#' @description
#' Creates shell commands to convert UMI counts from PRO-Seq data into positive strand and negative strand BigWig files.
#'
#' @param umi.counts Path to the input UMI counts file. Only a single file is allowed.
#' @param output.prefix Prefix for the output BigWig files. If not provided, it is derived from the input UMI counts filename.
#' @param bw.output.folder Directory for the BigWig files. Default= "db/bw/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file labels ("ps.bw", "ns.bw").
#' - `path`: Paths to the positive strand and negative strand BigWig files.
#' - `cmd`: Shell command to run the UMI to BigWig conversion.
#'
#' @examples
#' # Convert UMI counts to BigWig files for PRO-Seq data
#' cmd <- cmd_umiToBigwigProseq(
#'   umi.counts = "/data/umi_counts/sample_UMI_counts.txt"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_umiToBigwigProseq <- function(umi.counts,
                                  output.prefix= NULL,
                                  bw.output.folder= "db/bw/",
                                  Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if umi.counts file exists to allow wrapping!) ----
  if(length(umi.counts)!=1)
    stop("A unique umi.counts file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".txt$", "", basename(umi.counts))

  # Output file ----
  ps.bw <- file.path(bw.output.folder, paste0(output.prefix, ".ps.bw"))
  ns.bw <- file.path(bw.output.folder, paste0(output.prefix, ".ns.bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "umiToBigwigProseq.R", package = "vlite"),
    umi.counts,
    output.prefix
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("ps.bw", "ns.bw"),
                    path= c(ps.bw, ns.bw),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
