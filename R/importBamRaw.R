#' Import raw bam using samtools and fread
#'
#' Uses samtools to import a bam file.
#' @param file bam file path.
#' @param extra_arg Extra arg to be passed to samtools view.
#' @param headN Number of starting lines to import.
#' @param samtools_exe_path Path to samtools executable.
#'
#' @return Imported bam file
#' @export
#'
#' @examples
#' importBamRaw("path/to/bam/file.bam")
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
