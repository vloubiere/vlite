devtools::load_all()

set.seed(1)
a <- randomRegionsBSgenome("dm6", widths= sample(seq(250, 750), 10000, replace = T))
set.seed(2)
b <- randomRegionsBSgenome("dm6", widths= sample(seq(250, 750), 10000, replace = T))

t0 <- Sys.time()
closestBed(a, b)
t1 <- Sys.time()
t1-t0

a <- GRanges(a)
b <- GRanges(b)
test <- suppressWarnings(nearest(a, b, select= "all"))
a <- importBed(a[queryHits(test)])
b <- b[subjectHits(test)]
