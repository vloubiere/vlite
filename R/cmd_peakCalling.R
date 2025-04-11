#' Generate Commands for Peak Calling Using MACS2
#'
#' @description
#' Creates shell commands to perform peak calling on BAM files using MACS2.
#' Supports single-sample and merged peak calling, with optional input controls,
#' and outputs peak files and bedGraph files.
#'
#' @param bam A character vector specifying path(s) to treatment BAM file(s).
#' @param bam.input A character vector specifying path(s) to input BAM file(s). Default: `NULL`.
#' @param layout Sequencing layout, either `"SINGLE"` or `"PAIRED"`.
#' @param output.prefix Prefix for the output files.
#' @param keep.dup Number of duplicates to keep or `"all"`. Default: `1`.
#' @param extsize Integer value by which reads should be extended.
#'    Typically 200 for ChIP-Seq, DNAse, CutNrun; 75 for ATAC-Seq. Default: `200`.
#' @param shift Integer value by which reads will be shifted.
#'    Typically 0 for ChIP-Seq, DNAse & CutNrun; -35 for ATAC-Seq, -100 for DNAse. Default: `0`.
#' @param genome.macs2 Genome size parameter for MACS2 (e.g., `"mm"`, `"hs"`).
#' @param peaks.output.folder Directory for peak files. Default: `"db/peaks/"`.
#' @param broad Logical. Whether to call broad peaks. Default: `FALSE`.
#'
#' @return A `data.table` with:
#' - `file.type`: Output file labels (`"peaks"`, `"bedgraph"`).
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run the peak calling pipeline.
#'
#' @examples
#' # Basic peak calling without input controls
#' cmd <- cmd_peakCalling(
#'   bam = "/data/bam/rep1.bam",
#'   layout = "PAIRED",
#'   genome.macs2 = "mm",
#'   output.prefix = "sample1"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Peak calling with input controls
#' cmd <- cmd_peakCalling(
#'   bam = "/data/bam/rep1.bam",
#'   bam.input = "/data/bam/input1.bam",
#'   layout = "PAIRED",
#'   genome.macs2 = "mm",
#'   output.prefix = "sample1"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_peakCalling <- function(bam,
                            bam.input= NULL,
                            layout,
                            output.prefix,
                            keep.dup= 1,
                            extsize= 200,
                            shift= 0,
                            genome.macs2,
                            peaks.output.folder= "db/peaks/",
                            broad= FALSE)
{
  # Check ----
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")
  if(keep.dup!="all" && keep.dup %% 1!=0)
    stop("keep.dup should either be an integer value or set to 'all'")
  if(!is.na(extsize) && extsize %% 1!=0)
    stop("extsize should either be an integer value or NA (in which case, peak model will be computed).")
  if(shift %% 1!=0)
    stop("shift should be an integer.")
  if(!genome.macs2 %in% c("hs", "mm", "dm", "ce"))
    stop("genome.macs2 should be one of 'hs', 'mm', 'dm', 'ce'. See MACS2 documentation.")

  # Output files paths ----
  peaks.file <- file.path(peaks.output.folder,
                          paste0(output.prefix, ifelse(broad, "_peaks.broadPeak", "_peaks.narrowPeak")))
  bdg.file <- file.path(peaks.output.folder,
                        paste0(output.prefix, "_treat_pileup.bdg"))

  # Command ----
  cmd <- paste("macs2 callpeak -B --SPMR -g", genome.macs2, # Normalize tracks
               "-t", paste(bam, collapse = " "), # treatment bam file(s)
               "--keep-dup", keep.dup, # Keep duplicates
               "--outdir", peaks.output.folder, # Output directory
               "-n", output.prefix) # Output name
  if(layout=="PAIRED") # Format
    cmd <- paste(cmd, "-f BAMPE")
  if(!is.null(bam.input)) # Input bam file(s)
    cmd <- paste(cmd, "-c", paste(bam.input, collapse = " "))
  if(broad) # broad or narrow
    cmd <- paste(cmd, "--broad")
  if(!is.na(extsize)) # extsize (skip model)
    cmd <- paste(cmd, "--nomodel --extsize", extsize)
  if(shift!=0) # shift
    cmd <- paste(cmd, "--shift", shift)

  # Wrap commands output ----
  cmd <- data.table(file.type= c("peaks", "bedgraph"),
                    path= c(peaks.file, bdg.file),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
