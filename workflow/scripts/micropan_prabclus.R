# Pangenome openness and gene-content structure (Ortakci-group method):
# micropan::heaps (alpha < 1 => open pangenome) and prabclus::jaccard + classical PCoA (cmdscale),
# on the Roary gene_presence_absence.Rtab. micropan is not on conda; it is installed from CRAN on first use.
args <- commandArgs(trailingOnly = TRUE)
rtab <- args[1]; out <- args[2]; nperm <- as.integer(args[3])
if (!requireNamespace("micropan", quietly = TRUE)) {
  install.packages(c("microseq", "micropan"), repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressMessages({library(micropan); library(prabclus)})
set.seed(1)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
pa <- as.matrix(read.delim(rtab, row.names = 1, check.names = FALSE))
pan <- t(pa)                                   # genomes x gene families
h <- heaps(pan, n.perm = nperm)
write.table(data.frame(parameter = names(h), value = as.numeric(h)), file.path(out, "heaps.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
rf <- rarefaction(pan, n.perm = 1000)
write.table(rf, file.path(out, "rarefaction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
jd <- jaccard(pa); dimnames(jd) <- list(colnames(pa), colnames(pa))
pc <- cmdscale(as.dist(jd), k = 2, eig = TRUE)
var <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
write.table(data.frame(genome = rownames(pc$points), PC1 = pc$points[, 1], PC2 = pc$points[, 2], var1 = var[1], var2 = var[2]),
            file.path(out, "pcoa_jaccard.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.matrix(jd), file.path(out, "jaccard_distance.tsv"), sep = "\t", quote = FALSE)
cat(sprintf("Heaps alpha = %.3f (%s); PCoA %.1f%% / %.1f%%\n", h[2], ifelse(h[2] < 1, "open", "closed"), var[1], var[2]))
