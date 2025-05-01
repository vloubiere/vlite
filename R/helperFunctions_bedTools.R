# Import charcter strings coordinates
.importCharacterCoor <- function(coor)
{
  # Use the provided coordinates as names
  current <- data.table(name= coor)
  # Check if strand is provided
  stranded <- all(grepl("^.*:.*:.*$", coor))
  # Extract strand if relevant (and remove it from coordinates)
  if(stranded)
  {
    current[, strand:= gsub(".*:(.*)$", "\\1", name)]
    coor <- gsub("(.*):.*$", "\\1", coor)
  }
  # Remove big separator
  coor <- gsub(",", "", coor)
  # Coordinates data.table
  coorDT <- tstrsplit(coor, ":|-", type.convert = TRUE)
  coorDT <- as.data.table(coorDT)
  if(ncol(coorDT)==2)
  {
    # Add end if missing
    coorDT[, V3:= V2]
  }else if (ncol(coorDT)>3)
    stop("Character strings should follow the syntax: 'chr:start-end[:strand]'")
  
  # Name columns and cbind
  setnames(coorDT, c("seqnames", "start", "end"))
  current <- cbind(coorDT, current)
  
  # Return input type
  return(current)
}

# Import bed file
.importBedFile <- function(file, col.names= NULL)
{
  # Check format
  type <- unique(gsub(".*[.](.*)$", "\\1", file))
  if(length(type)>1)
    stop("Mixed input file formats are not supported.")
  
  # Import and rbind regions ----
  current <- lapply(file, data.table::fread)
  current <- rbindlist(current)
  
  # Column names ----
  if(is.null(col.names)) {
    col.names <- if(type=="bed") {
      c("seqnames", "start", "end", "name", "score", "strand")[1:ncol(current)]
    } else if(type=="broadPeak") {
      c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
    } else if(type=="narrowPeak") {
      c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
    }
  }
  setnames(current, col.names)
  
  # 0-base offset to 1-offset
  current[, start:= start+1]
  
  # Return
  return(current)
}

# Check formatting for bed file providd as data.table. See ?importBed()
.checkDataTableBedColClasses <- function(bed)
{
  # Use the provided coordinates as names
  current <- data.table::copy(bed)
  if(!all(c("seqnames", "start") %in% names(current)))
    stop("Genomic ranges should contain at least 'seqnames' and 'start' fields.")
  if(any(duplicated(names(current))))
     stop("Input column names should be unique.")

  # Check columns classes
  if ("seqnames" %in% names(current)) current[, seqnames := as.character(seqnames)]
  if ("start" %in% names(current)) current[, start := as.integer(start)]
  if ("end" %in% names(current)) {
    current[, end := as.integer(end)]
    # Check Boundaries
    if(any(current[, start>end]))
      warning("bed file contains genomic ranges with start>end -> malformed?")
  } 
  if ("name" %in% names(current)) current[, name := as.character(name)]
  if ("score" %in% names(current)) current[, score := as.numeric(score)]
  if ("strand" %in% names(current))
  {
    # Format
    current[, strand := as.character(strand)]
    current[strand==".", strand := "*"]
    # Check consistency
    if(any(!current$strand %in% c("+", "-", "*")))
      stop("All strand values should be one of c('+', '-' or '*') -> malformed genomic ranges!")
  }
  if ("signalValue" %in% names(current)) current[, signalValue := as.numeric(signalValue)]
  if ("pValue" %in% names(current)) current[, pValue := as.numeric(pValue)]
  if ("qValue" %in% names(current)) current[, qValue := as.numeric(qValue)]
  if ("peak" %in% names(current)) current[, peak := as.integer(peak)]

  # Add gr (genomic ranges) to class attributes
  setattr(x = current,
          name = "class",
          value = c("data.table", "data.frame", "gr"))

  # Return
  return(current)
}

# Bin a region using Nbins in ?binBed()
.bin_regions_using_Nbins <- function(regions,
                                     nbins)
{
  # Add nbins to regions ----
  regions[, nbins:= nbins]

  # Compute bins ----
  bins <- regions[, {

    # Order coordinates depending on strand ----
    coor <- if(neg.strand)
      end:start else
        start:end

    # Compute bins and indexes ----
    .c <- if(end-start+1>=nbins) {
      # Less bins than nucleotides
      bin.idx <- cut(seq_along(coor), breaks = nbins, labels = FALSE)
      .c <- data.table(coor, bin.idx)
      .c[, .(start= min(coor), end= max(coor)), bin.idx]
    } else {

      # More bins than nucleotides
      coor <- rep(coor, length.out= nbins)
      coor <- sort(coor, decreasing = neg.strand)
      bin.idx <- seq_along(coor)
      data.table(bin.idx, start= coor, end= coor)
    }

    # Ordering based on strand ----
    .c[order(bin.idx, decreasing = neg.strand)]
  }, line.idx]

  # Return ----
  return(bins)
}

# Bin a region using width parameters in ?binBed()
.bin_regions_using_width <- function(regions,
                                     bins.width,
                                     steps.width,
                                     bins.width.min)
{
  # Add bins.width and steps.width to regions ----
  regions[, bins.width:= bins.width]
  regions[, steps.width:= steps.width]

  # Compute bins ----
  bins <- regions[, {
    # Compute bins start and end depending on strand ----
    if(neg.strand){
      end <- seq(end, start, -steps.width)
      ext <- end-bins.width+1
      start <- ifelse(ext<start, start, ext)
    } else {
      start <- seq(start, end, steps.width)
      ext <- start+bins.width-1
      end <- ifelse(ext>end, end, ext)
    }

    # Make ranges DT ----
    .c <- data.table(bin.idx= seq_along(start),
                     start= as.integer(start),
                     end= as.integer(end))

    # Remove shorter bins ----
    if(bins.width.min)
      .c <- .c[end-start+1==bins.width]

    # Ordering based on strand ----
    .c[order(bin.idx, decreasing = neg.strand)]
  }, line.idx]

  # Return ----
  return(bins)
}
