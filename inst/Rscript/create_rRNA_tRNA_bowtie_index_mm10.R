ncRNA <- seqinr::read.fasta("/groups/stark/vloubiere/genomes/Mus_musculus/ncRNAs/Mus_musculus.GRCm38.ncrna.fa.gz")
annot <- sapply(ncRNA, attr, "Annot")
annot <- unlist(tstrsplit(annot, " ", keep= 5))
table(annot)

cmd <- c(
  "cd /groups/stark/vloubiere/genomes/Mus_musculus/ncRNAs/",
  "gunzip -c Mus_musculus.GRCm38.ncrna.fa.gz > Mus_musculus.GRCm38.ncrna.fa",
  "mkdir bowtie_index/",
  "bowtie-build Mus_musculus.GRCm38.ncrna.fa bowtie_index/mm10_ncrna"
)
cmd <- data.table(
  file.type= "bowtie_idx",
  path= "bowtie_index/mm10_ncrna",
  cmd= paste0(cmd, collapse = ";"),
  mem= 16,
  cores= 4,
  modules= list(c("build-env/2020", "bowtie/1.2.2-foss-2018b"))
)
vl_submit(cmd,
          logs = "/groups/stark/vloubiere/genomes/Mus_musculus/ncRNAs/logs/")
