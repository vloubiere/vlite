# Random sampling functions ----
rdmSamplingBS <- function(gSize,
                          prob,
                          widths)
{
  # Initiate data.table ----
  rdm <- data.table(seqnames= sample(x = gSize$seqnames,
                                     size = length(widths),
                                     prob = prob,
                                     replace = T),
                    width= widths)
  rdm[gSize, end:= i.end, on= "seqnames"]

  # Random sampling  ----
  rdm[, start:= sample(end, .N, replace = T), end]
  rdm[start>(end-width+1), start:= end-width+1]
  rdm[, end:= start+width-1]
  rdm[, strand:= "*"]

  # Check  ----
  if(any(rdm$start<1))
    stop("Error in rdmSamplingBS's helperFunction")

  # Return  ----
  return(rdm)
}
