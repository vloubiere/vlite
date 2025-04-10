# Download ATAC-Seq FC between PHD11 and control from the Nature paper
tmp <- tempfile(fileext = ".txt.gz")
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222193/suppl/GSE222193%5FATAC%5FTransient%5Fph%2DKD%5FD11%5Fvs%5FControl%5Fno%5Fph%2DKD%5FDESeq2%5FFC.txt.gz",
              destfile = tmp)

# Retrieve FC
FC <- fread(tmp)
FC <- cbind(importBed(FC$ID), FC[, !"ID"])

# Resize
res <- resizeBed(bed = FC,
                 center = "center",
                 upstream = 250,
                 downstream = 250,
                 genome = "dm6")

# Select JASPAR motifs
load("/groups/stark/vloubiere/motifs_db/vl_Dmel_motifs_DB_full.RData")
sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds]
pwms <- vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_perc]
mot.cluster <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_cluster]

# Compute motif counts
counts <- vl_motifCounts(bed= res,
                         genome= "dm6",
                         pwm_log_odds= sel)

# Lasso regression
lasso <- motifLassoRegression(response= FC$log2FoldChange,
                              counts = counts,
                              names= mot.cluster)

# Check PCC of the test set
vl_par()
plot(lasso$test_set)
addPcc(lasso$PCC)

# Retrieve top predictors
preds <- lasso$best_predictors

# Select top predictors
preds <- preds[, .SD[which.max(abs(s0))], name]
top <- preds[abs(s0)>0.03]
setorderv(top, "s0")

# Plot
barplot(top$s0,
        horiz= TRUE,
        names.arg = top$name,
        xlab= "Lasso coefficient")
