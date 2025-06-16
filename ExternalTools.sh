#!/bin/bash

# Create a directory to store tools
TOOLS_DIR="$HOME/rnaseq_tools"
mkdir -p "$TOOLS_DIR"
cd "$TOOLS_DIR"

echo "Downloading and installing required tools..."

# 1. sratoolkit
if ! command -v fastq-dump &> /dev/null; then
    echo "Downloading sratoolkit..."
    SRATOOLKIT_VERSION="3.0.7-ubuntu64"
    wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz -O sratoolkit.tar.gz
    tar -xzf sratoolkit.tar.gz
    SRATOOLKIT_DIR=$(find . -maxdepth 1 -type d -name "sratoolkit*" | head -n 1)
else
    echo "sratoolkit already installed."
fi

# 2. fastqc
if ! command -v fastqc &> /dev/null; then
    echo "Downloading FastQC..."
    wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -O fastqc.zip
    gunzip -q fastqc.zip
    chmod +x FastQC/fastqc
else
    echo "FastQC already installed."
fi

# 3. fastp
if ! command -v fastp &> /dev/null; then
    wget http://opengene.org/fastp/fastp
    chmod a+x ./fastp
else
    echo "fastp already installed."
fi

# 4. htseq
if ! command -v htseq-count &> /dev/null; then
    echo "Installing htseq via pip..."
    pip install --user htseq
else
    echo "htseq already installed."
fi

# 5. STAR
if ! command -v STAR &> /dev/null; then
    echo "Downloading STAR..."
    STAR_VERSION="2.7.10a"
    wget -q https://github.com/alexdobin/STAR/archive/refs/tags/$STAR_VERSION.tar.gz -O STAR.tar.gz
    tar -xzf STAR.tar.gz
    # Build STAR
    cd STAR-$STAR_VERSION/source
    make STAR
    cd ../../
else
    echo "STAR already installed."
fi

# Add binaries to PATH for this session
echo "Adding tools to PATH..."

## Add tools to PATH Variable
# TOOLS_DIR="$HOME/rnaseq_tools"
# SRATOOLKIT_DIR=$(find $TOOLS_DIR -maxdepth 1 -type d -name "sratoolkit*" | head -n 1)
# export PATH="$SRATOOLKIT_DIR/bin:$TOOLS_DIR/FastQC:$TOOLS_DIR:$TOOLS_DIR/STAR-2.7.10a/bin/Linux_x86_64_static:$PATH"
