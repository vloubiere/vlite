#
#
#
# covBam <- function(bed,
#                    bam,
#                    bedtools_path= "/software/2020/software/bedtools/2.27.1-foss-2018b/bin/bedtools") {
#   # Compose the bedtools command
#   cmd <- paste(bedtools_path, "coverage -a", bed, "-b", bam)
#
#   # Read the output directly into R
#   cov_dt <- data.table::fread(cmd = cmd, header = FALSE)
#
#   # The last column is the number of overlapping reads per region
#   return(cov_dt[[length(cov_dt)]])
# }
#
# bed.file <- "/groups/stark/vloubiere/projects/sebastian/db/bed/ORFtag/CpG1_input_rep1.1_mm10_collapsed_unique_insertions.bed"
# bed <- importBed(bed.file)
# bam.file <- "/groups/stark/vloubiere/projects/sebastian/db/bam/ORFtag/CpG1_input_rep1.1_mm10_collapsed.bam"
# bam <- importBam(bam.file,
#                  sel = c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq"))
# bed$cov <- covBam(bed= bed.file,
#                   bam = bam.file)
