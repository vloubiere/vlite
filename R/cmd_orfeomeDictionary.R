#' Generates ORFeome dictionary
#'
#' @description
#' Generates a shell command to identify confident ORFs/BCs combinations.
#'
#' @param BC.fqs A vector of paths to .fq.gz files containing the trimmed BC sequences.
#' @param bam.ORF Path to a .bam file containing reads aligned to ORF sequences.
#' @param minNreads Minimum number of supporting reads for a BC/ORF combination to be considered valid.
#' @param output.prefix Prefix for output files.
#' @param output.folder Output directory for result files. Default: "db/dictionary/".
#' @param Rpath Path to the Rscript binary. Default: "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file label ("peaks", "stats").
#' - `path`: Path to the output file.
#' - `cmd`: Shell command to run the peak calling pipeline.
#' - `job.name`: Default name for the job = "Dictionary".
#'
#' @examples
#' cmd <- cmd_orfeomeDictionary(
#'   BC.fqs = "/data/exp.fq.gz",
#'   bam.ORF = "/data/aligned.bam"
#' )
#' vl_submit(cmd, execute = FALSE)
#'
#' @export
cmd_orfeomeDictionary <- function(BC.fqs,
                                  bam.ORF,
                                  minNreads= 3,
                                  output.prefix,
                                  output.folder= "db/dictionary/",
                                  cores= 8,
                                  Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if files exist to allow wrapping!) ----
  if(!all(grepl(".fq.gz$|.fastq.gz$", BC.fqs)))
    stop("BC.fqs should be in .fq.gz format.")
  if(anyDuplicated(BC.fqs))
    stop("Some paths in BC.fqs are duplicated.")
  if(length(bam.ORF)>1 || !grepl(".bam$", bam.ORF))
    stop("bam.ORF should be a unique .bam file.")
  if(!is.numeric(minNreads) || minNreads %% 1 !=0)
    stop("minNreads should be an integer value.")

  # Output files paths ----
  stats.file <- file.path(output.folder, paste0(output.prefix, "_stats.txt"))
  stats.pdf <- file.path(output.folder, paste0(output.prefix, "_stats.pdf"))
  dictionary.rds <- file.path(output.folder, paste0(output.prefix, ".rds"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "ORFeome_make_dictionary.R", package = "vlite"),
    paste0(unique(BC.fqs), collapse= ","), # 1/ A comma-separated list of .fq.gz files containing the BC sequences
    bam.ORF, # 2/ A bam file containing ORF alignments for all the reads present in the BC file (arg[1])
    minNreads, # 3/ Minimum number of supporting reads
    file.path(output.folder, output.prefix) # 4/ Output prefix for output file
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("stats.file", "stats.pdf", "dictionary.rds"),
                    path= c(stats.file, stats.pdf, dictionary.rds),
                    cmd= cmd,
                    cores= cores,
                    job.name= "Dictionary")

  # Return ----
  return(cmd)
}
