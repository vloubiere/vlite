#' Convert Multiple BigWig Files to a Merged BigWig Format
#'
#' @description
#' Creates shell commands to merge multiple BigWig files and convert the merged BedGraph to BigWig format using a genome sizes file.
#'
#' @param bw Character vector of paths to input BigWig files. Only a single file is allowed; if more are provided, only the first is used and a warning is issued.
#' @param output.prefix Prefix for the output merged BigWig file.
#' @param bw.output.folder Directory for the output files. Default is "db/bw/".
#' @param genome A BSgenome ("mm10", "dm3"...).
#' @param Rpath Path to the Rscript binary. Default is "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A `data.table` with:
#' - `file.type`: Output file label ("bw").
#' - `path`: Path to the merged BigWig file.
#' - `cmd`: Shell command to run the BigWig merge and BedGraph to BigWig conversion.
#'
#' @examples
#' # Merge BigWig files and convert to a single BigWig
#' cmd <- cmd_bedgraphToBigwig(
#'   bw = c("/data/bw/sample_1.bw", "/data/bw/sample_2.bw"),
#'   output.prefix = "sample",
#'   genome = "mm10"
#' )
#' vl_submit(cmd, execute = FALSE)
#'
#' @export
cmd_mergeBigwig <- function(bw,
                            output.prefix,
                            bw.output.folder= "db/bw/",
                            genome,
                            Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks
  bw <- unique(bw)
  if(length(bw)==1)
    warning("Only one bigwig file was provided for the merge...")

  # Output file ----
  bdg.file <- file.path(bw.output.folder, paste0(output.prefix, "_merge.bdg"))
  bw.file <- file.path(bw.output.folder, paste0(output.prefix, "_merge.bw"))

  # Command ----
  cmd <- paste(
    "/software/2020/software/kent_tools/20190507-linux.x86_64/bin/bigWigMerge",
    paste0(bw, collapse = " "),
    bdg.file)
  cmd1 <- paste(
    Rpath,
    system.file("Rscript", "bedgraph_to_bigwig.R", package = "vlite"),
    bdg.file,
    bw.file,
    genome
  )
  cmd <- paste0(cmd, "; ", cmd1, "; rm ", bdg.file)

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= bw.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
