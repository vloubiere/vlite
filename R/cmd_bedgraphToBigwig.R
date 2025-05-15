#' Convert BedGraph to BigWig Format
#'
#' @description
#' Creates shell commands to convert a BedGraph file to BigWig format using a genome sizes file.
#'
#' @param bdg Path to the input BedGraph file. Only a single file is allowed.
#' @param output.prefix Prefix for the output BigWig file. If missing, will be constructed from the
#' basename of the bdg file.
#' @param bw.output.folder Directory for the BigWig file. Default= "db/bw/".
#' @param genome A BSgenome ("mm10", "dm3"...).
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A `data.table` with:
#' - `file.type`: Output file label ("bw").
#' - `path`: Path to the BigWig file.
#' - `cmd`: Shell command to run the BedGraph to BigWig conversion.
#'
#' @examples
#' # Convert a BedGraph file to BigWig format
#' cmd <- cmd_bedgraphToBigwig(
#'   bdg = "/data/peaks/sample.bdg",
#'   output.prefix = "sample",
#'   genome = "mm10"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_bedgraphToBigwig <- function(bdg,
                                 output.prefix,
                                 genome,
                                 scaling.factor= 1L,
                                 bw.output.folder= "db/bw/",
                                 Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if bdg file exists to allow wrapping!) ----
  if(length(bdg)!=1)
    stop("Only one bedgraph file should be provided.")
  if(missing(output.prefix))
    output.prefix <- gsub(".bdg$|.bedgraph$", "", basename(bdg))

  # Output file ----
  bw.file <- file.path(bw.output.folder, paste0(output.prefix, ".bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "bedgraph_to_bigwig.R", package = "vlite"),
    bdg, # Input bedgraph file
    bw.file, # Output bigwig file
    genome, # BSgenome used to construct chrom.sizes
    scaling.factor # Used to divide bedgraph scores (1= no scaling)
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= bw.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
