FROM python:3.11-slim

# System deps + R + headers needed to compile Seurat from source.
# Seurat compile is slow (~15 min); cache this layer aggressively.
RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base r-base-dev \
        build-essential \
        libcurl4-openssl-dev libxml2-dev libssl-dev \
        libhdf5-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Only Seurat + Matrix are needed for the RData -> MTX bridge
# (src/ingest/seurat_rdata.py).
RUN R -e "install.packages(c('Matrix','Seurat'), repos='https://cloud.r-project.org')"

WORKDIR /app

COPY pyproject.toml README.md ./
COPY wombat/ wombat/
COPY src/ src/
COPY scripts/ scripts/
COPY configs/ configs/

RUN pip install --no-cache-dir .

ENTRYPOINT ["wombat"]
