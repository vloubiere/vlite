#' Convert Multiple BigWig Files to a Merged BigWig Format
#'
#' @description
#' Creates shell commands to merge multiple BigWig files and convert the merged BedGraph to BigWig format using a
#' genome sizes file.
#'
#' @param bw Character vector of paths to input BigWig files.
#' @param output.prefix Prefix for the output merged BigWig file. If missing, will be constructed from the
#' basename of the first provided bw path.
#' @param bw.output.folder Directory for the output files. Default is "db/bw/".
#' @param tmp.bdg.folder Folder to save temporary .bedgraph files. Default is "db/bw/tmp/".
#' @param genome A BSgenome ("mm10", "dm3"...).
#' @param Rpath Path to the Rscript binary. Default is "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param bigWigMergePath Path to the tool executable.
#' Default= "/software/2020/software/kent_tools/20190507-linux.x86_64/bin/bigWigMerge"
#'
#' @return A `data.table` with:
#' - `file.type`: Output file label ("bw").
#' - `path`: Path to the merged BigWig file.
#' - `cmd`: Shell command to run the BigWig merge and BedGraph to BigWig conversion.
#' - `job.name`: Default name for the job = "bwMerge".
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
                            tmp.bdg.folder= "db/bw/tmp/",
                            genome,
                            Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                            bigWigMergePath= "/software/2020/software/kent_tools/20190507-linux.x86_64/bin/bigWigMerge")
{
  # Checks
  bw <- unique(bw)
  if(length(bw)==1)
    warning("Only one bigwig file was provided and will be copied as is.")
  if(missing(output.prefix))
    output.prefix <- gsub(".bw$", "", basename(bw[1]))

  # Output file ----
  bdg.file <- file.path(tmp.bdg.folder, paste0(output.prefix, "_merge.bdg"))
  bw.file <- file.path(bw.output.folder, paste0(output.prefix, "_merge.bw"))

  # Command to merge bigwigs to bedgraph ----
  merge.cmd <- paste("mkdir", dirname(bdg.file), ";",
                     bigWigMergePath, # Path to the tool
                     paste0(bw, collapse = " "),
                     bdg.file)

  # Transforme merged bedgraph into a bigwig ----
  cmd <- cmd_bedgraphToBigwig(bdg = bdg.file,
                              output.prefix = gsub(".bw$", "", basename(bw.file)),
                              bw.output.folder = bw.output.folder,
                              genome = genome,
                              Rpath = Rpath)

  # If several files were provided ----
  if(length(bw)>1) {
    # Add merging command
    cmd[, cmd:= paste0(merge.cmd, "; ", cmd)]
    # Remove bdg.file
    cmd[, cmd:= paste0(cmd, "; rm ", bdg.file)]
  } else {
    # Simply copy the input file
    cmd$cmd <- paste("cp", bw, bw.file)
  }

  # Return ----
  cmd[, job.name= "bwMerge"]
  return(cmd)
}
