#' Call Peaks from bigWig Files
#'
#' @description
#' Generates a shell command to call peaks from a ChIP-seq experiment bigWig file, with optional input control, over specified genomic regions.
#' Outputs a .narrowPeak file and a statistics .txt file.
#'
#' @param experiment.bw.file Path to the experiment bigWig (.bw) file. Must be unique.
#' @param input.bw.file Optional path to input/control bigWig (.bw) file. Must be unique or NULL.
#' @param bed.subset Optional path to a bed file. If specified, peak calling will be restricted to these regions.
#' Note that the strand is not considered. Default= NULL.
#' @param output.prefix Prefix for output files. If missing, it will be inferred from the name of the experiment.bw.file.
#' @param output.folder Output directory for result files. Default: "db/peaks/".
#' @param zscore.cutoff Z-score cutoff for putative peaks. Default: 1.64.
#' @param local.enr.cutoff Local background enrichment cutoff. Default: 1.5.
#' @param local.fdr.cutoff Local background FDR cutoff. Default: 0.05.
#' @param input.enr.cutoff Input enrichment cutoff. Default: 1.5.
#' @param input.fdr.cutoff Input FDR cutoff. Default: 0.05.
#' @param Rpath Path to the Rscript binary. Default: "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file label ("peaks", "stats").
#' - `path`: Path to the output file.
#' - `cmd`: Shell command to run the peak calling pipeline.
#' - `job.name`: Default name for the job = "peakCalling".
#'
#' @examples
#' cmd <- cmd_peakCallingFromBw(
#'   experiment.bw.file = "/data/exp.bw",
#'   input.bw.file = "/data/input.bw"
#' )
#' vl_submit(cmd, execute = FALSE)
#'
#' @export
cmd_peakCallingFromBw <- function(experiment.bw.file,
                                  input.bw.file= NULL,
                                  bed.subset= NULL,
                                  output.prefix,
                                  output.folder= "db/peaks/",
                                  zscore.cutoff= 1.64,
                                  local.enr.cutoff= 1.5,
                                  local.fdr.cutoff= 0.05,
                                  input.enr.cutoff= 1.5,
                                  input.fdr.cutoff= 0.05,
                                  Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if files exist to allow wrapping!) ----
  if(length(experiment.bw.file)!=1)
    stop("A unique experiment bw file should be provided.")
  if(length(input.bw.file)>1)
    stop("If provided, the input bw file should be unique.")
  if(!all(grepl(".bw$", c(experiment.bw.file, input.bw.file))))
    stop("Experiment and optional input bw file should be in .bw extension.")
  if(!is.null(bed.subset) && (length(bed.subset)>1 || !grepl(".bed$", bed.subset)))
    stop("If specified, bed should be a path to a unique .bed file.")
  if(missing(output.prefix))
    output.prefix <- gsub(".bw$", "_peaks", basename(experiment.bw.file))
  if(!is.numeric(c(zscore.cutoff, local.enr.cutoff, local.fdr.cutoff, input.enr.cutoff, input.fdr.cutoff)))
    stop("zscore.cutoff, local.enr.cutoff, local.fdr.cutoff, input.enr.cutoff and input.fdr.cutoff should be numeric.")


  # Output files paths ----
  peaks.file <- file.path(output.folder, paste0(output.prefix, ".narrowPeak"))
  stats.file <- file.path(output.folder, paste0(output.prefix, "_stats.txt"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "peakCallingFromBw.R", package = "vlite"),
    experiment.bw.file, # Experiment bw file
    ifelse(is.null(input.bw.file), "NULL", input.bw.file), # Optional input bw file
    ifelse(is.null(bed.subset), "NULL", bed.subset), # An optional bed file restricting peak calling to certain regions
    peaks.file, # Output .narrowPeak file path
    stats.file, # Output .txt statistics file
    zscore.cutoff, # z-score cutoff to identify putative peaks. Default= 1.64
    local.enr.cutoff, # Local background enrichment cutoff. Default= 1.5
    local.fdr.cutoff, # Local background FDR cutoff. Default= 0.05
    input.enr.cutoff, # Input enrichment cutoff. Default= 1.5
    input.fdr.cutoff # Input FDR cutoff. Default= 0.05
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("peaks", "stats"),
                    path= c(peaks.file, stats.file),
                    cmd= cmd,
                    job.name= "peakCalling")

  # Return ----
  return(cmd)
}
