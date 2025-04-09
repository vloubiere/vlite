# Download EpiCancer clusters
tmp <- tempfile(pattern = ".xlsx")
download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07328-w/MediaObjects/41586_2024_7328_MOESM5_ESM.xlsx",
              destfile = tmp)

# Import genes
genes <- readxl::read_xlsx(tmp, sheet = 1)
genes <- as.data.table(genes)

# Gene clusters
cl <- genes[!cluster %in% c("NA", "Unaffected")]
cl <- split(cl$FBgn, cl[, .(cluster, V2= ifelse(PcG_bound, "Bound", "Unbound"))])

# Compare activators and repressor Cellular Compartments
enr <- vl_GOenrich(geneIDs = cl,
                   geneUniverse.IDs = genes$FBgn,
                   species= "Dm")
enr[, cl:= factor(cl,
                  c("Reversible.Bound",
                    "Irreversible.Bound",
                    "Transient-specific.Bound",
                    "Down 1.Bound",
                    "Down 2.Bound",
                    "Down 3.Bound",
                    "Reversible.Unbound",
                    "Irreversible.Unbound",
                    "Transient-specific.Unbound",
                    "Down 1.Unbound",
                    "Down 2.Unbound",
                    "Down 3.Unbound"))]

# Plot
vl_par(mai= c(.9, 2, .9, 1.3))
plot(obj= enr,
     top.enrich = 5,
     padj.cutoff= 0.05,
     order = "log2OR",
     cex= .5)
