#' Generate Commands for Peak Calling Using MACS2
#'
#' @description
#' Generates shell commands to perform peak calling on BAM files using MACS2,
#' converts bedGraph to bigWig format, and identifies confident peaks across replicates.
#' The function handles both single-sample and merged peak calling, with optional input controls.
#'
#' @param bam character. Vector of paths to BAM files for peak calling.
#' @param bam.input character. Optional vector of paths to input/control BAM files.
#'        Must match the order of treatment BAM files. Default: NULL.
#' @param layout character(1). Sequencing layout, must be either "PAIRED" or "SINGLE".
#' @param output.prefix character(1). Prefix for output files.
#' @param keep.dup
#' @param extsize Default= 200.
#' @param shift
#' @param genome.macs2 character(1). Genome size parameter for MACS2 (e.g., "mm" or "hs").
#' @param peaks.output.folder character(1). Directory where peak files will be written.
#' @param Rpath character(1). Path to Rscript executable.
#'        Default: "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"
#' @param genome.bw character(1). Genome sizes file for bigWig conversion.
#' @param bw_output_folder character(1). Directory where bigWig files will be written.
#' @param extsize numeric(1). Read extension size for MACS2. Default: 300.
#' @param broad logical(1). Whether to call broad peaks. Default: FALSE.
#'
#' @return A data.table with three columns:
#' \describe{
#'   \item{file.type}{Labels for output files ("peaks", "bw", "peaks.merge", "bw.merge", "conf.peaks")}
#'   \item{path}{Full paths to the output files}
#'   \item{cmd}{Shell command to run the peak calling pipeline}
#' }
#'
#' @details
#' The function generates a pipeline that:
#' \enumerate{
#'   \item Calls peaks for individual replicates using MACS2
#'   \item Calls peaks on merged BAM files
#'   \item Converts MACS2 bedGraph output to bigWig format
#'   \item Identifies confident peaks across replicates
#' }
#'
#' @section Output Files:
#' The function generates several types of files:
#' \itemize{
#'   \item Individual replicate peaks: <prefix>_<replicate>_peaks.narrowPeak
#'   \item Individual replicate signal: <prefix>_<replicate>.bw
#'   \item Merged peaks: <prefix>_merge_peaks.narrowPeak
#'   \item Merged signal: <prefix>_merge.bw
#'   \item Confident peaks: <prefix>_confident_peaks.narrowPeak
#' }
#'
#' @section Requirements:
#' \itemize{
#'   \item MACS2 must be installed and available in the system PATH
#'   \item R with required packages (genomicsPipelines) must be installed
#'   \item Output directories must exist and be writable
#' }
#'
#' @examples
#' \dontrun{
#' # Basic peak calling without input controls
#' cmd_peakCallling(
#'   bam = c("/data/bam/rep1.bam", "/data/bam/rep2.bam"),
#'   replicate = c("rep1", "rep2"),
#'   layout = "PAIRED",
#'   genome.macs2 = "mm",
#'   output.prefix = "sample1",
#'   peaks.output.folder = "/data/peaks/",
#'   genome.bw = "/data/genome/mm10.chrom.sizes",
#'   bw_output_folder = "/data/bigwig/"
#' )
#'
#' # Peak calling with input controls
#' cmd_peakCallling(
#'   bam = c("/data/bam/rep1.bam", "/data/bam/rep2.bam"),
#'   replicate = c("rep1", "rep2"),
#'   bam.input = c("/data/bam/input1.bam", "/data/bam/input2.bam"),
#'   layout = "PAIRED",
#'   genome.macs2 = "mm",
#'   output.prefix = "sample1",
#'   peaks.output.folder = "/data/peaks/",
#'   genome.bw = "/data/genome/mm10.chrom.sizes",
#'   bw_output_folder = "/data/bigwig/"
#' )
#' }
#'
#' @seealso
#' \itemize{
#'   \item MACS2 documentation: \url{https://github.com/macs3-project/MACS}
#' }
#'
#' @export
cmd_peakCallling <- function(bam,
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
  if(length(bam)>1)
    stop("If multiple bam files are provided, their paths should be concatenated and space-separated.")
  if(!is.null(bam.input) && length(bam.input)>1)
    stop("If multiple bam.input files are provided, their paths should be concatenated and space-separated.")
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")
  if(keep.dup!="all" && keep.dup %% 1!=0)
    stop("keep.dup should either be an integer value or set to 'all'")
  if(!is.na(extsize) && extsize %% 1!=0)
    stop("extsize should either be an integer value or NA (in which case, peak model will be computed).")
  if(shift %% 1!=0)
    stop("shift should either be an integer.")
  if(!genome.macs2 %in% c("hs", "mm", "dm", "ce"))
    stop("genome.macs2 should be one of 'hs', 'mm', 'dm', 'ce'. See MACS2 documentation.")
  if(!dir.exists(peaks.output.folder))
    dir.create(peaks.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output files paths ----
  peaks.file <- file.path(peaks.output.folder,
                          paste0(output.prefix, ifelse(broad, "_peaks.broadPeak", "_peaks.narrowPeak")))
  bdg.file <- file.path(peaks.output.folder,
                        paste0(output.prefix, "_treat_pileup.bdg"))

  # Command ----
  cmd <- paste("macs2 callpeak -B --SPMR -g", genome.macs2, # Normalize tracks
               "-t", bam, # bam input file(s)
               "--keep-dup", keep.dup, # Keep duplicates
               "--outdir", peaks.output.folder, # Output directory
               "-n", output.prefix) # Output name
  if(layout=="PAIRED") # Format
    cmd <- paste(cmd, "-f BAMPE")
  if(!is.null(bam.input)) # Input bam file(s)
    cmd <- paste(cmd, "-c", bam.input)
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
