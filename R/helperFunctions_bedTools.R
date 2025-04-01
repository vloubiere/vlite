# Check formatting for bed file providd as data.table. See ?importBed()
.checkDataTableBedColClasses <- function(bed)
{
  # Use the provided coordinates as names
  current <- data.table::copy(bed)

  # Check columns classes
  if ("seqnames" %in% names(current)) current[, seqnames := as.character(seqnames)]
  if ("start" %in% names(current)) current[, start := as.integer(start)]
  if ("end" %in% names(current)) current[, end := as.integer(end)]
  if ("name" %in% names(current)) current[, name := as.character(name)]
  if ("score" %in% names(current)) current[, score := as.numeric(score)]
  if ("strand" %in% names(current))
  {
    current[, strand := as.character(strand)]
    current[strand==".", strand := "*"]
    # Check condisten strand
    if(any(!current$strand %in% c("+", "-", "*")))
      stop("All strand values should be one of c('+', '-' or '*') -> malformed ranges!")
  }

  # Check Boundaries
  if("end" %in% names(current) && any(current[, start>end], na.rm = T))
    warning("bed file contains ranges with start>end -> malformed ranges!")

  # Return
  return(current)
}

# Import bed file
.importBedFile <- function(file, extra.columns= NULL)
{
  # Import and rbind regions
  current <- lapply(file, data.table::fread)
  current <- rbindlist(current)

  # Add names
  col.names <- c("seqnames", "start", "end", "name", "score", "strand", extra.columns)[1:ncol(current)]
  setnames(current, col.names)

  # 0-base offset to 1-offset
  current[, start:= start+1]

  # Return
  return(current)
}

# Import narrowPeak file. See ?importBed()
.importNarrowPeakFile <- function(file)
{
  # Import and rbind regions
  current <- lapply(file,
                    data.table::fread,
                    col.names= c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"),
                    colClasses= c("character", "integer", "integer", "character", "numeric", "character", "numeric", "numeric", "numeric", "integer"))
  current <- rbindlist(current)

  # 0-base offset to 1-offset
  current[, start:= start+1]

  # Return
  return(current)
}

# Import broadPeak file. See ?importBed()
.importBroadPeakFile <- function(file)
{
  # Import and rbind regions
  current <- lapply(file,
                    data.table::fread,
                    col.names= c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue"),
                    colClasses= c("character", "integer", "integer", "character", "numeric", "character", "numeric", "numeric", "numeric"))
  current <- rbindlist(current)

  # 0-base offset to 1-offset
  current[, start:= start+1]

  # Return
  return(current)
}

# Check input type for ?importBed()
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
  {
    stop("Format could not be determined.
           If you are using character coordinates,
           use a consistent syntax for all regions: 'chrX:1000-2000' or 'chrX:1000-2000'")
  }

  # Name columns and cbind
  setnames(coorDT, c("seqnames", "start", "end"))
  current <- cbind(coorDT, current)

  # Return input type
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
