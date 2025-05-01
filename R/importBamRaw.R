#' Import Raw BAM File Using Samtools and fread
#'
#' Uses samtools to import a BAM file and reads the output into R using fread.
#'
#' @param file A character string specifying the path to the BAM file.
#' @param extra_arg Optional extra arguments that will be passed to `samtools view`.
#' @param headN An integer specifying the number of starting lines to import. By default, all reads will be processed.
#' @param samtools_exe_path A character string specifying the path to the samtools executable.
#' Default= "/software/2020/software/samtools/1.9-foss-2018b/bin/samtools".
#'
#' @return
#' A data table containing the imported BAM file.
#'
#' @examples
#' # Import a BAM file
#' importBamRaw("path/to/bam/file.bam")
#'
#' # Import with additional arguments
#' importBamRaw("path/to/bam/file.bam", extra_arg = "-f 4", headN = 100)
#'
#' @export
importBamRaw <- function(file,
                         extra_arg,
                         headN,
                         samtools_exe_path= "/software/2020/software/samtools/1.9-foss-2018b/bin/samtools")
{
  # Initiate command
  cmd <- paste(samtools_exe_path, "view")

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
