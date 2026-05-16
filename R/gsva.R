#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA is not installed; cannot compute GSVA/ssGSEA scores")
  }
  if (!requireNamespace("GSEABase", quietly = TRUE)) {
    stop("GSEABase is not installed; cannot read GMT files")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript R/gsva.R <expression_csv> <gmt_file> <output_csv>")
}

expression_path <- args[[1]]
gmt_path <- args[[2]]
output_path <- args[[3]]

expr <- read.csv(expression_path, check.names = FALSE, row.names = 1)
gene_sets <- GSEABase::getGmt(gmt_path)

if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  params <- GSVA::ssgseaParam(as.matrix(expr), gene_sets)
  scores <- GSVA::gsva(params)
} else {
  scores <- GSVA::gsva(as.matrix(expr), gene_sets, method = "ssgsea", verbose = FALSE)
}

write.csv(scores, output_path)
