#' Convert BedGraph to BigWig Format
#'
#' @description
#' Creates shell commands to convert a BedGraph file to BigWig format using a genome sizes file.
#'
#' @param bdg Path to the input BedGraph file. Only a single file is allowed.
#' @param output.prefix Prefix for the output BigWig file.
#' @param bw.output.folder Directory for the BigWig file. Default: `"db/bw/"`.
#' @param genome Path to the genome sizes file.
#' @param Rpath Path to the Rscript binary. Default: `"/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"`.
#'
#' @return A `data.table` with:
#' - `file.type`: Output file label (`"bw"`).
#' - `path`: Path to the BigWig file.
#' - `cmd`: Shell command to run the BedGraph to BigWig conversion.
#'
#' @examples
#' # Convert a BedGraph file to BigWig format
#' cmd <- cmd_bedgraphToBigwig(
#'   bdg = "/data/peaks/sample.bdg",
#'   output.prefix = "sample",
#'   genome = "/data/genome/mm10.chrom.sizes"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_bedgraphToBigwig <- function(bdg,
                                 output.prefix,
                                 bw.output.folder= "db/bw/",
                                 genome,
                                 Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(length(bdg)!=1)
    stop("Only one bedgraph file should be provided.")
  if(!dir.exists(bw.output.folder))
    dir.create(bw.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output file ----
  bw.file <- file.path(bw.output.folder, paste0(output.prefix, "_", genome, ".bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "bedgraph_to_bigwig.R", package = "vlite"),
    bdg,
    bw.file,
    genome
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= bw.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
