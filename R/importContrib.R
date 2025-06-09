#' Import contributions object
#'
#' @param h5 A vector of paths to h5 file(s) containing the contribution scores.
#' @param bed Bed file(s) containing the regions for which contributions were computed. By default, the file is searched in the same folder as the h5 file.
#' @param fa fasta file(s) containing the sequences for which contributions were computed. By default, the file is searched in the same folder as the h5 file.
#' @param FUN The function to apply to each layer of the h5 object to retrieve the contribution scores. By default, each layer is a
#' matrix of dimensions nchar(sequence)*4. Thus, default= function(x) rowSums(x).
#'
#'
#' @return A contribution data.table containing, for each region, the contribution scores (as list) and the corresponding sequence.
#' @export
#'
#' @examples
#' folder <- "/groups/stark/vloubiere/projects/epiDeepCancer/db/model_PH18/contributions/"
#' files <- list.files(folder, "h5$", recursive = TRUE, full.names = TRUE)
#' dat <- importContrib(h5= files)
#'
importContrib <- function(h5,
                          bed= list.files(dirname(h5), ".bed$", full.names = TRUE),
                          fa= list.files(dirname(h5), ".fa$", full.names = TRUE),
                          FUN= function(x) rowSums(x))
{
  # Metadata ----
  meta <- data.table(h5= h5,
                     bed= bed,
                     fa= fa)

  # Import ----
  dat <- meta[, {
    # bed regions
    .c <- importBed(bed)
    if(!is.na(fa) && file.exists(fa))
      .c$seq <- as.character(seqinr::read.fasta(fa, as.string = TRUE)) else
        message("No valid fasta file found")
    if(!is.na(h5) && file.exists(h5))
    {
      .h <- rhdf5::h5read(h5, "contrib_scores/class")
      dim <- paste0(dim(.h), collapse = "*")
      print(paste("h5 dimensions:", dim))
      .c$score <- lapply(seq(dim(.h)[3]), function(i) FUN(.h[,,i]))
    }else
      message("No valid h5 file found")
    .c
  }, .(h5, bed, fa)]
  dat$h5 <- dat$bed <- dat$fa <- NULL

  # Message ----
  if(any(nchar(dat$seq)>20))
  {
    options(datatable.prettyprint.char = 10)
    message("Printing option set to 10. To reset it to default, use options(datatable.prettyprint.char = NULL)")
  }

  # Return object ----
  return(dat)
}
