# bulkRNAseq_python
A Python-based pipeline for downloading, processing, and analyzing bulk RNA sequencing (RNA-seq) data from the Gene Expression Omnibus (GEO). This toolkit automates key steps including data retrieval, quality control, alignment, counting, normalization, and downstream analysis to facilitate reproducible RNA-seq workflows.

## Features
Data Download: Automates retrieval of raw RNA-seq data from GEO using fastq_dump.py.

Quality Control: Performs QC checks on FASTQ files with fastq_qc.py.

Reference Genome Preparation: Builds reference indices for alignment (mkref.py).

Alignment: Uses STAR aligner (STAR.py) to map reads to the reference genome.

Counting: Quantifies gene expression levels using htseq.py.

Normalization: Normalizes raw counts to prepare data for analysis (Normalize.py).

Downstream Analysis: Supports differential expression and visualization workflows (DownstreamAnalysis.py).

Utility Tools: Helper functions for file handling, process management, and logging (OS_Tools.py, RunTools.py).

## Installation
Clone the repository:

```
git clone https://github.com/ewalsh3963/bulkRNAseq_python.git
cd bulkRNAseq_python
```

Install dependencies (recommended via conda or pip):
```
pip install -r requirements.txt
```

Install external command line tools reuqired for pipeline (assumes Ubuntu system)

- fastq_dump (sratoolkit)
- fastqc
- fastp
- STAR
- htseq-count

```
cd bulkRNAseq_python
chmod +x ExternalTools.sh
./ExternalTools.sh

## Add Variables to PATH
TOOLS_DIR="$HOME/rnaseq_tools"
SRATOOLKIT_DIR=$(find $TOOLS_DIR -maxdepth 1 -type d -name "sratoolkit*" | head -n 1)
export PATH="$SRATOOLKIT_DIR/bin:$TOOLS_DIR/FastQC:$TOOLS_DIR:$TOOLS_DIR/STAR-2.7.10a/bin/Linux_x86_64_static:$PATH"
```

## Usage
Run individual pipeline steps as needed. Example commands:

```
python fastq_dump.py --accession GSEXXXXX --outdir ./data
python fastq_qc.py --input ./data
python mkref.py --genome fasta/genome.fa --outdir ./ref
python STAR.py --read1 ./data/sample_1.fastq --read2 ./data/sample_2.fastq --outdir ./aligned
python htseq.py --bam ./aligned/sample.bam --gtf annotations.gtf --out ./counts
python Normalize.py --counts ./counts --outdir ./normalized
python DownstreamAnalysis.py --input ./normalized --outdir ./results
## (Replace arguments with your specific data paths and options.)
```
## Logging
The pipeline logs key actions and parameters for reproducibility. Logging utilities are implemented to create and append detailed logs for each step.

## Contributing
Contributions and improvements are welcome! Please open issues or pull requests on GitHub.

