#' Find Motif Positions in Sequences or Genomic Regions
#'
#' @description
#' A wrapper around motifmatchr::matchMotifs that maps motif positions in a set of genomic coordinates or sequences or using a PWMatrixList as input.
#'
#' @param sequences A named character vector of sequences to analyze. This argument takes precedence over the bed argument.
#' @param bed Genomic ranges in a format compatible with ?importBed, from which genomic sequences will be retrieved when sequences is set to NULL.
#' @param pwm_log_odds A PWMatrixList (in log2 odds ratio format) containing motifs to map. For example, see "/groups/stark/vloubiere/motifs_db/".
#' @param genome The genome to use as background when bg = "genome" and/or to retrieve sequences (when bed is specified).
#' Default= NULL.
#' @param bg The background model for motif detection. Options are "genome", "subject" (inferred from input sequences) or "even" (0.25, 0.25, 0.25, 0.25). Default= "genome".
#' @param p.cutoff The p-value cutoff for motif detection. Default= 5e-5.
#' @param pos.strand If set to TRUE, only motifs on the positive strand are considered (default= FALSE).
#' @param collapse.overlapping Logical. If TRUE, motifs with 70 percent overlap ore more will be collapsed, irrespective of their strand. Default= FALSE.
#' @param scratch The path to the scratch directory for temporary files. Default= "/scratch-cbe/users/vincent.loubiere/motifs/".
#' @param cleanup.cache If set to TRUE, overwrites cached intermediate results. Default= FALSE.
#'
#' @return A data.table containing motif positions for each sequence, with the following columns:
#' - `motif`: The name of the motif.
#' - `seqlvls`: The sequence names.
#' - `mot.count`: The number of motifs detected.
#' - `ir`: A nested `data.table` with columns for motif positions (`start`, `end`, `strand`, `width`, `score`).
#'
#' @examples
#' # Load Drosophila motifs
#' pfm.file <- system.file("extdata/hand_curated_Dmel_motifs_SCENIC_lite_Dec_2025.pfm", package = "vlite")
#' pwm <- importJASPAR(pfm.file)$pwms_log_odds
#'
#' # Test sequence containing two overalapping Jra motifs plus a non overlapping one
#' test <- c("Jra"= "CATGAGTCAGGTGAGAAATGAGTCAT")
#' mot <- vl_motifPos(sequences = test, pwm_log_odds = pwm["Jra"], genome = "dm3", p.cutoff = 5e-4)
#' motifPosToBed(mot)
#' coll <- vl_motifPos(sequences = test, pwm_log_odds = pwm["Jra"], genome = "dm3", p.cutoff = 5e-4, collapse.overlapping = T)
#' motifPosToBed(coll)
#'
#' # Download Dev and Hk enhancers from pe-STARR-Seq paper
#' tmp <- tempfile(pattern = '.xlsx')
#' download.file(
#' url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-52921-2/MediaObjects/41467_2024_52921_MOESM4_ESM.xlsx',
#' destfile = tmp
#' )
#' enh <- as.data.table(readxl::read_xlsx(tmp, sheet = 2, skip = 1))
#' set <- rbind(
#'   enh[group=='dev' & detail %in% c('medium', 'strong')][1:100, .(group, seqnames, start, end, enhancer_sequence)],
#'   enh[group=='hk' & detail %in% c('medium', 'strong'), .(group, seqnames, start, end, enhancer_sequence)]
#' )
#' seq <- set$enhancer_sequence
#' names(seq) <- make.unique(set$group)
#'
#' # Map motifs
#' mot <- vl_motifPos(sequences = seq, pwm_log_odds = pwm[c("Dref", "Jra")], bg = "subject", p.cutoff = 5e-4)
#'
#' # Look at counts: Jra should be more frequent in developemtan and Dref in housekeeping
#' mot[, class:= tstrsplit(seqlvls, "[.]", keep= 1)]
#' mot[, rm.NA:= ifelse(is.na(mot.count), 0, mot.count)]
#' vl_boxplot(rm.NA~motif+class, mot, outline= T)
#'
#' # Convert to scores matrix (useful for heatmaps)
#' mat <- vlite::motifPosToMatrix(mot, seqWidth = nchar(seq)[1])$Dref
#' vl_par()
#' vl_heatmap(
#'   mat,
#'   cluster.rows = set$group,
#'   show.row.clusters = "left",
#'   main= "Dref max scores",
#'   show.colnames = FALSE,
#'   legend.title = "Score"
#' )
#'
#' @export
vl_motifPos <- function(sequences,
                        bed,
                        pwm_log_odds,
                        genome= NULL,
                        bg= "genome",
                        p.cutoff= 5e-5,
                        pos.strand= FALSE,
                        collapse.overlapping= FALSE,
                        scratch= "/scratch-cbe/users/vincent.loubiere/motifs/",
                        cleanup.cache= FALSE)
{
  # Checks ----
  if(!missing(sequences) && !missing(bed))
    warning("sequences are provided -> input bed will not be used.")
  if(missing(sequences) && missing(bed))
    stop("sequences of bed regions should be specified.")
  if(is.null(genome) && (missing(sequences) | bg=="genome"))
    stop("genome is set to NULL.")
  if(!"PWMatrixList" %in% class(pwm_log_odds))
    pwm_log_odds <- do.call(TFBSTools::PWMatrixList, pwm_log_odds)
  if(anyDuplicated(sapply(pwm_log_odds, TFBSTools::name)))
    stop("Duplicated motif names in provided PWMs. Check them using TFBSTools::name().")
  if(!is.numeric(p.cutoff) || p.cutoff>1)
    stop("p.cutoff should be a numeric value <= 1")

  # Get sequences ----
  if(missing(sequences)) {
    bed <- importBed(bed)
    sequences <- getBSsequence(bed, genome)
  }

  # Make sure sequence names are unique ----
  if(is.null(names(sequences)) || anyDuplicated(names(sequences)))
    stop("All sequences should have a unique name!")

  # Create tmp output folder ----
  params <- list(sequences,
                 ifelse(bg=="genome", genome, bg),
                 bg,
                 p.cutoff)
  sub.folder <- digest::digest(params)
  tmp.folder <- paste0(scratch, "/", sub.folder, "/")
  print(paste0("Temp files will be stored in '", tmp.folder, "'"))
  dir.create(tmp.folder,
             showWarnings = FALSE,
             recursive = TRUE)

  # Create final file cache name ----
  final.file <- vl_cache_file(
    list(
      sub.folder,
      pwm_log_odds,
      pos.strand,
      collapse.overlapping
    ),
    tmp.dir = tmp.folder,
    extension = ".rds"
  )

  # Main function ----
  if(!file.exists(final.file) | cleanup.cache) {

    # Retrieve unique motif names ----
    mot.names <- sapply(pwm_log_odds, TFBSTools::name)
    mot.names <- make.unique(mot.names)

    # Generate output file names ----
    output.files <- paste0(tmp.folder, "/", gsub("/", "__", mot.names), ".rds")

    # Check for existing files ----
    missing.files <- !file.exists(output.files) | cleanup.cache

    # Map motifs ----
    if(any(missing.files))
    {
      print(paste0(sum(!missing.files), "/", length(pwm_log_odds), " motif files already existed!"))

      # Call position ----
      parallel::mcmapply(function(mot, output.file)
      {
        res <- motifmatchr::matchMotifs(pwms = mot,
                                        subject = sequences,
                                        genome = genome,
                                        p.cutoff = p.cutoff,
                                        bg = bg,
                                        out = "positions")[[1]]
        saveRDS(res, output.file)
      },
      mot= pwm_log_odds[missing.files],
      output.file= output.files[missing.files],
      mc.preschedule = TRUE,
      mc.cores = max(c(1, data.table::getDTthreads()-1)))
    }
    print("All motif positions computed ;)")

    # Post-processing ----
    pos <- data.table(motif= mot.names,
                      file= output.files)
    final <- pos[, {
      # Import
      .c <- readRDS(file)
      .c <- as.data.table(.c)

      # Select positive strand motifs
      if(pos.strand)
        .c <- .c[strand=="+"]

      # Add seqnames
      .c[, seqnames:= names(sequences)[group]]

      # If specified, collapsed motifs that overlap >70%, ignore.strand
      if(collapse.overlapping && nrow(.c))
      {
        # Collapse
        coll <-  collapseBed(bed = .c, ignore.strand = TRUE)

        # Re-split large regions into bins corresponding to 70% motif width
        mot.width <- ceiling(mean(.c$width)*0.7)
        coll$nBins <- floor(coll[, end-start+1]/mot.width)
        if(any(coll$nBins>1)) {
          coll[, coll.idx:= .I]
          coll <- rbind(
            coll[nBins==1],
            coll[nBins>1, binBed(.SD, nBins), nBins],
            fill= T
          )
          setorderv(coll, c("coll.idx", "bin.idx"))
          coll$coll.idx <- coll$line.idx <- coll$bin.idx <- NULL
        }
        coll$nBins <- NULL

        # Compute width
        coll[, width:= end-start+1]

        # Retrieve score
        coll$score <- .c[coll, max(score), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1

        # Overwrite current table
        .c <- coll
      }

      # Simplify ir to list
      .c[, seqlvls:= seqnames]
      cols <- intersect(names(.c), c("start", "end", "strand", "width", "score"))
      .c <- .c[, .(mot.count= .N, ir= .(.SD)), seqlvls, .SDcols= cols]

      # Add missing levels
      all <- data.table(seqlvls= names(sequences))
      res <- merge(all, .c, by= "seqlvls", all.x= TRUE, sort= FALSE)

      # Return
      res
    }, motif]

    # Save cache file
    saveRDS(final,
            final.file)
  } else {
    message("All intermediate files exist -> importing final table.")
    final <- readRDS(final.file)
  }


  # Return ----
  return(final)
}
