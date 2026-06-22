#!/usr/bin/env Rscript
#
# Install the R dependencies used by the R helper scripts in this directory.
# Targeted against R 4.4 (release) with Bioconductor 3.19+; genefu, GSVA, and
# GSEABase are Bioconductor packages, mice is on CRAN.
#
# Required packages, and the script that uses each:
#   genefu    -- R/pam50.R         (official PAM50 subtyping + rorS / ROR-S)
#   GSVA      -- R/gsva.R          (ssGSEA pathway scoring)
#   GSEABase  -- R/gsva.R          (getGmt: read .gmt gene-set files)
#   mice      -- R/mice_impute.R   (multiple imputation by chained equations)
#
# Usage:  Rscript R/install.R

options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- CRAN packages -----------------------------------------------------------
cran_pkgs <- c("BiocManager", "mice")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing CRAN package: %s", pkg))
    install.packages(pkg)
  }
}

# --- Bioconductor packages ---------------------------------------------------
bioc_pkgs <- c("genefu", "GSVA", "GSEABase")
to_install <- bioc_pkgs[!vapply(
  bioc_pkgs,
  function(p) requireNamespace(p, quietly = TRUE),
  logical(1)
)]
if (length(to_install) > 0) {
  message(sprintf("Installing Bioconductor packages: %s",
                  paste(to_install, collapse = ", ")))
  BiocManager::install(to_install, ask = FALSE, update = FALSE)
}

# --- Verify ------------------------------------------------------------------
all_pkgs <- c(cran_pkgs, bioc_pkgs)
missing <- all_pkgs[!vapply(
  all_pkgs,
  function(p) requireNamespace(p, quietly = TRUE),
  logical(1)
)]
if (length(missing) > 0) {
  stop(sprintf("Failed to install: %s", paste(missing, collapse = ", ")))
}
message("All R dependencies installed: ", paste(all_pkgs, collapse = ", "))
