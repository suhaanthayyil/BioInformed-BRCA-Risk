#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("genefu", quietly = TRUE)) {
    stop("genefu is not installed; cannot compute official PAM50/ROR scores")
  }
  library(genefu)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript R/pam50.R <expression_samples_x_genes_csv> <annotation_csv> <output_csv> [do_mapping]")
}

expression_path <- args[[1]]
annotation_path <- args[[2]]
output_path <- args[[3]]
do_mapping <- ifelse(length(args) >= 4, tolower(args[[4]]) %in% c("true", "1", "yes"), FALSE)

data(pam50, package = "genefu")
data(pam50.robust, package = "genefu")

expr <- read.csv(expression_path, check.names = FALSE, row.names = 1)
annot <- read.csv(annotation_path, check.names = FALSE)
if ("probe" %in% colnames(annot)) {
  rownames(annot) <- annot$probe
}

common <- intersect(colnames(expr), rownames(annot))
if (length(common) < 40) {
  stop(sprintf("Only %d PAM50 genes available; refusing official PAM50 scoring", length(common)))
}
expr <- expr[, common, drop = FALSE]
annot <- annot[common, , drop = FALSE]

sub <- molecular.subtyping(
  sbt.model = "pam50",
  data = as.matrix(expr),
  annot = annot,
  do.mapping = do_mapping,
  verbose = FALSE
)
ror <- rorS(
  data = as.matrix(expr),
  annot = annot,
  do.mapping = do_mapping,
  verbose = FALSE
)

out <- data.frame(
  sample_id = rownames(expr),
  pam50_subtype_genefu = as.character(sub$subtype),
  rorS_score = as.numeric(ror$score),
  rorS_risk = as.character(ror$risk),
  n_pam50_genes_input = ncol(expr),
  stringsAsFactors = FALSE
)
write.csv(out, output_path, row.names = FALSE)
