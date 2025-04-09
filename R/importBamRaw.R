#' Import Raw BAM File Using Samtools and fread
#'
#' Uses `samtools` to import a BAM file and reads the output into R using `fread`.
#'
#' @param file A character string specifying the path to the BAM file.
#' @param extra_arg Additional arguments to pass to `samtools view`. Optional.
#' @param headN An integer specifying the number of starting lines to import. Optional.
#' @param samtools_exe_path A character string specifying the path to the `samtools` executable.
#' Defaults to `/software/2020/software/samtools/1.9-foss-2018b/bin/samtools view`.
#'
#' @return
#' A data table containing the imported BAM file.
#'
#' @examples
#' \dontrun{
#' # Import a BAM file
#' importBamRaw("path/to/bam/file.bam")
#'
#' # Import with additional arguments
#' importBamRaw("path/to/bam/file.bam", extra_arg = "-f 4")
#'
#' # Import only the first 100 lines
#' importBamRaw("path/to/bam/file.bam", headN = 100)
#' }
#'
#' @export
importBamRaw <- function(file,
                         extra_arg,
                         headN,
                         samtools_exe_path= "/software/2020/software/samtools/1.9-foss-2018b/bin/samtools view")
{
  # Initiate command
  cmd <- samtools_exe_path

  # Add extra arguments
  if(!missing(extra_arg))
    cmd <- paste(cmd,
                 extra_arg)
  cmd <- paste(cmd, file)

  # Add head if specified
  if(!missing(headN))
  {
    if(!is.integer(headN))
      headN <- as.integer(headN)
    cmd <- paste(cmd, "| head -n", headN)
  }

  # Import
  fread(cmd= cmd,
        fill= T,
        header= F,
        sep= "\t")
}
