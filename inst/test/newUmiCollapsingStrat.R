# Find duplicated UMIs (for which another UMI with 1 mismatch has been sequenced >= 10 times more) ----
reads <- data.table(coor= "A",
                    UMI= c("GGGGGGGGGG", "CGGGGGGGGG", "TGGGGGGGGG", "GGGGGGCGGG", "CGGGGGCGGG", "CTGGGGGGGG"),
                    umi_N= c(1000, 500, 40, 100, 9, 40),
                    total_count= 1000)
cols <- paste0("dup.", seq(10))
for(i in 1:10) {
  reads[, (cols[i]):= umi_N<=umi_N[1]/10, .(coor, gsub(paste0("^(.{", i-1, "})."), "\\1", UMI))]
}
reads <- reads[rowSums(reads[, dup.1:dup.10])==0, !c(cols), with= FALSE]
reads
