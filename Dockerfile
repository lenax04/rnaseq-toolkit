# rnaseq-toolkit Docker image
# Provides a fully reproducible environment for RNA-seq analysis.
#
# Build:
#   docker build -t rnaseq-toolkit:0.1.0 .
#
# Run:
#   docker run --rm -v $(pwd)/data:/workspace/data \
#              -v $(pwd)/results:/workspace/results \
#              rnaseq-toolkit:0.1.0 \
#              rnaseq-toolkit --counts data/counts.csv \
#                             --metadata data/metadata.csv \
#                             --design "~condition" \
#                             --output results/

FROM python:3.12-slim-bookworm

LABEL maintainer="RNAseqKit Contributors <contact@rnaseqkit.org>"
LABEL version="0.1.0"
LABEL description="Streamlined RNA-seq differential expression and pathway enrichment"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libhdf5-dev \
    libopenblas-dev \
    liblapack-dev \
    libffi-dev \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /workspace

# Copy package files
COPY pyproject.toml setup.py README.md ./
COPY src/ ./src/

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir \
        pandas>=1.5.0 \
        numpy>=1.23.0 \
        pydeseq2>=0.4.1 \
        matplotlib>=3.7.1 \
        seaborn>=0.12.2 \
        gseapy>=1.0.5 \
        scikit-learn>=1.2.0 \
        statsmodels>=0.14.0 \
        scipy>=1.10.0 \
        GEOparse>=2.0.0 \
        snakemake>=7.25.0 && \
    pip install --no-cache-dir -e .

# Create workspace directories
RUN mkdir -p /workspace/data /workspace/results /workspace/logs

# Default entrypoint
ENTRYPOINT ["rnaseq-toolkit"]
CMD ["--help"]
