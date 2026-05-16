#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("genefu", quietly = TRUE)) {
    stop("genefu is not installed; cannot compute official PAM50/ROR scores")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript R/pam50.R <expression_csv> <annotation_csv> <output_csv>")
}

expression_path <- args[[1]]
annotation_path <- args[[2]]
output_path <- args[[3]]

expr <- read.csv(expression_path, check.names = FALSE, row.names = 1)
annot <- read.csv(annotation_path, check.names = FALSE)

res <- genefu::molecular.subtyping(
  sbt.model = "pam50",
  data = as.matrix(expr),
  annot = annot,
  do.mapping = TRUE
)

write.csv(res, output_path, row.names = FALSE)
