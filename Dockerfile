# Reference reproduction environment for BRCA-PathwayML (BMC AI revision, E3).
#
# Provides Python 3.12 (pinned deps from requirements.txt) + system R with the
# genefu/GSVA/mice packages needed for the official PAM50-ROR scoring path.
#
# NOTE: this image is heavy and slow to build -- compiling the Bioconductor
# `genefu` stack (and its dependencies) from source under r-base takes a while.
# It is intended as a documented, fully-specified environment for reviewers, not
# as a lightweight runtime. The Tier-1 reproduction (`make reproduce-headline`)
# does not actually need R; R is only required for the genefu PAM50-ROR
# comparator and the all-from-raw rebuild.
#
# Build:  docker build -t brca-pathwayml .
# Run:    docker run --rm -it brca-pathwayml

FROM python:3.12-slim

# System R + build/runtime deps for the Bioconductor stack and scientific wheels.
RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base \
        build-essential \
        gfortran \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        zlib1g-dev \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Python dependencies first (cached unless requirements.txt changes).
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# R dependencies (genefu/GSVA/GSEABase via Bioconductor, mice via CRAN).
COPY R/install.R ./R/install.R
RUN Rscript R/install.R

# Project sources.
COPY . .

CMD ["bash", "-lc", "echo 'BRCA-PathwayML reference environment.' && \
echo 'Reproduce the headline tables + TNBC endpoint with:' && \
echo '    make reproduce-headline' && \
echo 'Run the test suite (includes genefu PAM50-ROR) with:' && \
echo '    make test'"]
