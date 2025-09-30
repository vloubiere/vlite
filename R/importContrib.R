#' Import contributions object
#'
#' @param h5 A vector of paths to h5 file(s) containing the contribution scores.
#' @param bed Optional path(s) to Bed file(s) containing the coordinates for which contributions were computed.
#' @param fa Optional path(s) to fasta file(s) containing the sequences for which contributions were computed.
#' @param FUN The function to apply to each layer of the h5 object to retrieve the contribution scores.
#' By default, each layer is a matrix of dimensions nchar(sequence)*4. Thus, default= function(x) rowSums(x).
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
                          bed= NULL,
                          fa= NULL,
                          FUN= function(x) rowSums(x))
{
  # Checks ----
  if(any(!file.exists(h5)))
    stop("Some h5 file(s) were not found.")
  if(!is.null(bed) && any(!file.exists(bed)))
    stop("Some bed file(s) not found.")
  if(any(!file.exists(fa)) || length(fa)==0)
    stop("Some fasta file(s) not found.")

  # Metadata ----
  meta <- data.table(h5= h5,
                     bed= if(is.null(bed)) NA else bed,
                     fa= if(is.null(fa)) NA else fa)

  # Import ----

  dat <- meta[, {

    # Import h5
    .c <- rhdf5::h5read(h5, "contrib_scores/class")
    dim <- paste0(dim(.c), collapse = "*")
    print(paste("h5 dimensions:", dim))
    .c <- data.table(contrib.score= lapply(seq(dim(.c)[3]), function(i) FUN(.c[,,i])))

    # Add bed coordinates if provided
    if(!is.na(bed))
      .c <- cbind(.c, importBed(bed))

    # Add fasta sequences if provided
    if(!is.na(fa)) {
      # Import
      .fa <- seqinr::read.fasta(fa, as.string = TRUE)
      if(is.na(bed)) {
        # Create coordinates based on fasta
        .c[, seqnames:= names(.fa)]
        .c[, start:= 1]
        .c[, end:= nchar(as.character(.fa))]
      } else {
        # Fasta sequence name
        .c[, fa.name:= names(.fa)]
        if("name" %in% names(.c) && !all.equal(.c$name, .c$fa.name))
          warning("bed file and fasta names are different.")
      }
      # DNA sequence
      .c[, seq:= as.character(.fa)]
      if(any(nchar(.c$seq)>20)) {
        options(datatable.prettyprint.char = 50)
        message("Printing option set to 10. To reset it to default, use options(datatable.prettyprint.char = NULL)")
      }
    }

    # Add Dummy coordinates if necessary
    if(is.na(bed) && is.na(fa)) {
      .c[, seqnames:= paste0("seq", .I)]
      .c[, start:= 1]
      .c[, end:= sapply(contrib.score, length)]
    }

    # Return
    .c
  }, .(h5, bed, fa)]

  # Set Put score at the end, clean and save ----
  setcolorder(dat, "contrib.score", after = ncol(dat))
  dat$h5 <- dat$bed <- dat$fa <- NULL
  return(dat)
}
