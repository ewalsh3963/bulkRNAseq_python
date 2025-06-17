# !/bin/python3
# Download fastq files from Gene Expression Omnibus and perform RNA-seq analysis
# QBDS
# Evan Walsh (evanwalsh396@gmail.com)

import argparse
import logging
import os
from datetime import datetime
import re
import sys
import time
import glob
import mkref
import htseq
import Normalize
import STAR
import fastq_dump
import fastq_qc
import pandas as pd
import json
import subprocess
import DownstreamAnalysis
import pdb
import OS_Tools


class CommandLine:
    """ Handle the command line, usage and help requests """

    def __init__(self, in_opts=None):
        """ CommandLine constructor: Implements a parser to interpret the command line argv string using argparse. """
        self.parser = argparse.ArgumentParser(description="%(prog)s is the pipeline for downloading RNAseq data with SRA Toolkit. \n"
        "Modules in this pipeline include fastqDump, mkref, STAR, HTseq i\n",
        epilog='CAPSIDA BIOTHERAPEUTICS STAR alignment pipeline', 
        add_help=True, 
        prefix_chars='-', 
        usage='python3 %(prog)s [-h] [-v] [fullAlign, fastqDump, mkref, STAR, htseq, Normalzie, DownstreamAnalysis] ...')
     
        ## version 
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0 - Capsida Biotherapeutics - Evan Walsh (evanwalsh396@gmail.com)')

        subparsers = self.parser.add_subparsers(title="programs", 
                                                dest="program",
                                                metavar="[fullAlign, fastqDump, mkref, STAR, htseq, Normalzie]")
        ## Full Pipeline
        self.fullAlign_parse = subparsers.add_parser("fullAlign",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'fullAlign\' runs all modules to process fastq files into read counts:\n"
                                                )
        self.fullAlign_args()

        ## fastqDump
        self.fastqDump_parse = subparsers.add_parser("fastqDump",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'fastqDump\' download fastq files from Gene Expression Omnibus using SRA IDs")
        self.fastqDump_args()

        ## fastqDump
        self.fastqQC_parse = subparsers.add_parser("fastqQC",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'fastqQC\' run fastqc and adapter trimming on fastq files")
        self.fastqQC_args()

        ## mkref 
        self.mkref_parse = subparsers.add_parser("mkref",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'mkref\' is used to generate a STAR index for alignment of fastq files")
        self.mkref_args()

        ## STAR 
        self.STAR_parse = subparsers.add_parser("STAR",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'STAR\' alignment")
        self.STAR_args()
        
        ## htseq
        self.htseq_parse = subparsers.add_parser("htseq",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'htseq\' count reads from bam file")
        self.htseq_args()

        ## Normalize RNAseq counts
        self.Normalize_parse = subparsers.add_parser("Normalize",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'Normalize\' convert count data to normalized values")
        self.Normalize_args()

        ## Downstream Analysis
        self.DownstreamAnalysis_parse = subparsers.add_parser("DownstreamAnalysis",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'DownstreamAnalysis\' post-quantifcation analysis and normalization")
        self.DownstreamAnalysis_args()

        if in_opts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(in_opts)
                
    def fullAlign_args(self):
        self.fullAlign_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.fullAlign_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='/home/ewalsh/scratch')

        ## fastqDump parameters
        self.fullAlign_parse.add_argument('-srr', metavar='--srr_run_table', action='store', help='SRR Run table file with SRR IDs to use for download (full path)', required=False, type=str)
        self.fullAlign_parse.add_argument('-N', metavar='--N', action='store', help='Minimum spot id for read downsampling', required=False, default=None, type=str)
        self.fullAlign_parse.add_argument('-X', metavar='--X', action='store', help='Maximum spot id for read downsampling', required=False, default=None, type=str)

        ## STAR parameters
        self.fullAlign_parse.add_argument('-r', metavar='--reference_dir', action='store', required=False, type=str, default=None, help='path to the STAR index generated by split-pipe mkref')

        ## htseq parameters
        self.fullAlign_parse.add_argument("-htseq_args", type=json.loads, default={}, help="JSON string of additional htseq-count arguments")
        self.fullAlign_parse.add_argument('-gtf', metavar='--gtf_file', action='store', help='gtf file', required=False, type=str, default=None)
        self.fullAlign_parse.add_argument('-bed', metavar='--bed_file', action='store', help='gtf file', required=False, type=str, default=None)

    def fastqDump_args(self):
        self.fastqDump_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.fastqDump_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='/home/ewalsh/scratch')
        self.fastqDump_parse.add_argument('-srr', metavar='--srr_run_table', action='store', help='SRR Run table file with SRR IDs to use for download (full path)', required=True, type=str)

        ## Downsample parameters 
        self.fastqDump_parse.add_argument('-N', metavar='--N', action='store', help='Minimum spot id for read downsampling', required=False, default=None, type=str)
        self.fastqDump_parse.add_argument('-X', metavar='--X', action='store', help='Maximum spot id for read downsampling', required=False, default=None, type=str)

    def fastqQC_args(self):
        self.fastqQC_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.fastqQC_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='/home/ewalsh/scratch')
        self.fastqQC_parse.add_argument('-i', metavar='--input_dir', action='store', help='output file path', required=False, type=str, default=None)
  
    def mkref_args(self):
        self.mkref_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.mkref_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='/home/ewalsh/scratch')
        self.mkref_parse.add_argument('-species', metavar='--species', action='store', required=False, type=str, help='Species of samples [default: None]', default=None)
        
        self.mkref_parse.add_argument("-star_args", type=json.loads, default={}, help="JSON string of additional htseq-count arguments")

        ## reference files 
        self.mkref_parse.add_argument('-rf', metavar='--reference_fasta', action='store', required=False, type=str, help='genome reference fasta file', default=None)
        self.mkref_parse.add_argument('-rg', metavar='--reference_gtf', action='store', required=False, type=str, help='genome reference gtf file', default=None)

    def STAR_args(self):
        self.STAR_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.STAR_parse.add_argument('-si', metavar='--sr_Id', action='store', nargs = '*', type=str, required=False, help='SR IDs to align fastqs to')
        self.STAR_parse.add_argument('-r', metavar='--reference_dir', action='store', required=False, type=str, default=None, help='path to the STAR index generated by split-pipe mkref')
        
        self.STAR_parse.add_argument("-star_args", type=json.loads, default={}, help="JSON string of additional htseq-count arguments")

        self.STAR_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='/home/ewalsh/scratch')
        self.STAR_parse.add_argument('-i', metavar='--input_dir', action='store', help='output file path', required=False, type=str, default=None)

    def htseq_args(self):
        self.htseq_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.htseq_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='/home/ewalsh/scratch')
        self.htseq_parse.add_argument('-i', metavar='--input_dir', action='store', help='output file path', required=False, type=str, default=None)
        
        self.htseq_parse.add_argument("-htseq_args", type=json.loads, default={}, help="JSON string of additional htseq-count arguments")
        self.htseq_parse.add_argument('-gtf', metavar='--gtf_file', action='store', help='gtf file (must be unzipped)', required=False, type=str, default=None)
        self.htseq_parse.add_argument('-bed', metavar='--bed_file', action='store', help='bed file (must be unzipped)', required=False, type=str, default=None)

    def DownstreamAnalysis_args(self):
        self.DownstreamAnalysis_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.DownstreamAnalysis_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='/home/ewalsh/scratch')
        self.DownstreamAnalysis_parse.add_argument('-i', metavar='--input_dir', action='store', help='input file path', required=False, type=str, default=None)
        
        ## Filtering and curating data
        self.DownstreamAnalysis_parse.add_argument('-ms', metavar='--min_samples', action='store', help='minimum sample number for gene filtering', required=False, type=int, default=2)
        self.DownstreamAnalysis_parse.add_argument('-m', metavar='--metadata_file', action='store', help='sample meta data', required=False, type=str, default = None)
        
        ## DEseq2
        self.DownstreamAnalysis_parse.add_argument('-df', metavar='--design_factors', action='store', help='design factors', nargs = "*", required=False, type=str, default = None)
        self.DownstreamAnalysis_parse.add_argument('-c', '--contrast',  metavar=('factor', 'test', 'reference'),  nargs=3,  help='Contrast for DESeq2: format is factor test reference (e.g., Condition RS C)',  required=False,  type=str,  default=None)

        ## Stats
        self.DownstreamAnalysis_parse.add_argument('--padj_threshold', action='store', help='adj. p-value for significance threshold ', required=False, type=float, default = 0.05)
        self.DownstreamAnalysis_parse.add_argument('--abs_logFC', action='store', help='abs(log FC) for significance threshold ', required=False, type=float, default = 0.5)

    def Normalize_args(self):
        self.Normalize_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.Normalize_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='/home/ewalsh/scratch')
        self.Normalize_parse.add_argument('-i', metavar='--input_dir', action='store', help='input file path', required=False, type=str, default=None)
        
        ## Normalization parameters
        self.Normalize_parse.add_argument('-gtf', metavar='--gtf_file', action='store', help='gtf file', required=False, type=str, default=None)
        self.Normalize_parse.add_argument('-f', metavar='--format', action='store',help='Normalization format: one of RPKM, CPM, TPM, FPKM (default: None --> performs all operations)',required=False,type=str,default=None,choices=[None, 'RPKM', 'CPM', 'TPM', 'FPKM'])
        self.Normalize_parse.add_argument('-species', metavar='--species', action='store', required=False, type=str, help='Species of samples [default: None]', default=None)

class LogInitializer:
    """ Create Logs depending on the program chosen. """
    def __init__(self, program, out):
        self.program = program
        self.out = out

    def logger_check(self):
        """ Create and check logs. """
        run_logs = os.path.join(self.out, "{}_scanpy.log".format(self.program))  # create logs file name
        OS_Tools.check_file(run_logs, critical=False, rewrite=True)  # check logs file - create new log
        return run_logs

    def param_logger(self, short, run_logs, title, args):
        """ Add input parameters to logs. """
        breaker = '=' * 120
        my_logs = OS_Tools.Logs(run_logs)
        my_logs.append_logs("{} {} ran with the following parameters:\n{}\n{}\n{}\n\n"
                            .format("RunSpipe.py", self.program, breaker, "\n".join(f'{k} = {v}' for k, v in vars(args).items()),breaker))

        my_logs.append_logs("{}:\n{} | {}\n".format(short, title, datetime.now()))

def main(command=None):
    """ Run specified program(s). """
    if command is None:
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main

    ## check if outroot exists 
    my_command.args.o = os.path.abspath(my_command.args.o)
    OS_Tools.ensure_directory(my_command.args.o, critical=False) 
    outroot = my_command.args.o

    # Handle top-level options
    if my_command.args.program == "fullAlign":
        """ Run all steps of the RNA-seq pipeline """
        """ EXCEPT mkref... NEEDS TO BE RUN INDIVIDUALLY AND IS EXPECTED PRIOR TO RUNNING all COMMAND """

        init_log = LogInitializer("fastqDump", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
        init_log.param_logger("Running all alignment commands [fastq_dump, fastq_qc, STAR, htseq]", run_logs, my_command.args.s, my_command.args)

        ###################################################################################################
        ## 1) Run fastq-dump if a SRA run table is supplied
        my_logs.append_logs(breaker)
        my_logs.append_logs("Downloading fastq files from GEO using fastq_dump")
        outdir = os.path.join(outroot, 'expdata', my_command.args.s)
        ## Run fastq_dump.py
        my_command.args.root = outroot
        my_command.args.o = outdir
        fastq_dump.main(my_command.args, run_logs)

        ###################################################################################################
        ## 2) Run fastqc and adaptor trimming
        my_logs.append_logs(breaker)
        my_logs.append_logs("Running FastqQC and adapter trimming on fastq files")
        fastq_qc.main(my_command.args, run_logs)

        ###################################################################################################
        ## 3) Run STAR alignment
        my_logs.append_logs(breaker)
        my_logs.append_logs("Downloading fastq files from GEO using fastq_dump")

        outdir = os.path.join(outroot, 'results', my_command.args.s)
        my_command.args.o = outdir
        STAR.main(my_command.args, run_logs)

        ###################################################################################################
        ## 4) Run htseq_count
        my_logs.append_logs(breaker)
        my_logs.append_logs("Downloading fastq files from GEO using fastq_dump")

        outdir = os.path.join(outroot, 'PublicRNAseqData', my_command.args.s)
        my_command.args.o = outdir
        htseq.main(my_command.args, run_logs)

    elif my_command.args.program == "fastqDump":
        """ Download fastq files from GEO with fastq dump from SRA Toolkit."""
        outdir = os.path.join(outroot, 'expdata', my_command.args.s); OS_Tools.ensure_directory(outdir, critical=False) 
        init_log = LogInitializer("fastqDump", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
                
        init_log.param_logger("Download fastq files from GEO using fastq_dump", run_logs, my_command.args.s, my_command.args)
        
        ## Run fastq_dump.py
        my_command.args.root = my_command.args.o
        my_command.args.o = outdir
        fastq_dump.main(my_command.args, run_logs)

        ########### clean up output directory ###########
        cmd = ['rm', '-rf', "*_sra"]
        subprocess.run(cmd, cwd=outdir)

    elif my_command.args.program == "fastqQC":
        """ Generate a STAR reference for alignment of fastq files."""
        outdir = os.path.join(outroot, 'expdata', my_command.args.s); OS_Tools.ensure_directory(outdir, critical=False) 
        init_log = LogInitializer("fastqQC", outdir)
        run_logs = init_log.logger_check()

        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
        init_log.param_logger("Running FastqQC and adapter trimming on fastq files", run_logs, my_command.args.s, my_command.args)
        
        ## Run fastq_qc.py
        my_command.args.root = my_command.args.o
        my_command.args.o = outdir
        fastq_qc.main(my_command.args, run_logs)

    elif my_command.args.program == "mkref":
        """ Generate a STAR reference for alignment of fastq files."""

        ###################################################################################################
        ## 1) Make a STAR refrence for cellranger
        outdir = os.path.join(outroot, 'genomes', my_command.args.s); OS_Tools.ensure_directory(outdir, critical=False) 
        init_log = LogInitializer("mkref", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
                
        init_log.param_logger("Generate a STAR reference with fasta and gtf files", run_logs, my_command.args.s, my_command.args)
        ## Run mkref.py 
        my_command.args.o = outdir
        mkref.main(my_command.args, run_logs)

    elif my_command.args.program == "STAR":
        """ Run STAR alignment."""

        ###################################################################################################
        ## 1) Make a STAR refrence for cellranger
        outdir = os.path.join(outroot, 'results', my_command.args.s); OS_Tools.ensure_directory(my_command.args.o, critical=False) 
        init_log = LogInitializer("align", outdir)
        run_logs = init_log.logger_check()

        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
                
        init_log.param_logger("Run STAR alignment", run_logs, my_command.args.s, my_command.args)
        my_command.args.root = my_command.args.o
        my_command.args.o = outdir
        STAR.main(my_command.args, run_logs)

    elif my_command.args.program == "htseq":
        """ Run htseq-count"""

        ###################################################################################################
        ## 1) Make a STAR refrence for cellranger
        outdir = os.path.join(outroot, 'PublicRNAseqData', my_command.args.s); OS_Tools.ensure_directory(outdir, critical=False) 
        init_log = LogInitializer("htseq", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
                
        init_log.param_logger("htseq-count", run_logs, my_command.args.s, my_command.args)
        ## Run mkref.py 
        my_command.args.root = my_command.args.o
        my_command.args.o = outdir
        htseq.main(my_command.args, run_logs)
    
    elif my_command.args.program == "DownstreamAnalysis":
        """ Run QC, filtering, and normalization on post-quantifcation analysis"""

        ###################################################################################################
        ## 1) Make a STAR refrence for cellranger
        outdir = os.path.join(outroot, 'PublicRNAseqData', my_command.args.s); OS_Tools.ensure_directory(outdir, critical=False) 
        init_log = LogInitializer("DownstreamAnalysis", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
                
        init_log.param_logger("DownstreamAnalysis", run_logs, my_command.args.s, my_command.args)
        my_command.args.root = my_command.args.o
        my_command.args.o = outdir
        DownstreamAnalysis.main(my_command.args, run_logs)

    elif my_command.args.program == "Normalize":
        """ Run convery counts to RPKM values"""

        ###################################################################################################
        ## 1) Make a STAR refrence for cellranger
        outdir = os.path.join(outroot, 'PublicRNAseqData', my_command.args.s); OS_Tools.ensure_directory(outdir, critical=False) 
        init_log = LogInitializer("htseq", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
                
        init_log.param_logger("Normalize", run_logs, my_command.args.s, my_command.args)
        ## Run mkref.py 
        my_command.args.root = my_command.args.o
        my_command.args.o = outdir
        Normalize.main(my_command.args, run_logs)


if __name__ == "__main__":
    main()


# !import code; code.interact(local=vars())

###########################################################################
## CHO-1 GSE138292

# python3 /home/ewalsh/FastqProcess/main.py fastqDump -s CHO_1_GSE138292 \
# -srr /home/ewalsh/FastqProcess/RunFiles/CHO_1_GSE138292/SraRunTable.csv 

# python3 /home/ewalsh/FastqProcess/main.py fastqDump -s CHO_1_GSE138292 \
# -srr /home/ewalsh/FastqProcess/RunFiles/CHO_1_GSE138292/SraRunTable.csv  -N 1 -X 750000

# python3 /home/ewalsh/FastqProcess/main.py fastqQC -s CHO_1_GSE138292 

# python3 /home/ewalsh/FastqProcess/main.py mkref -s CHO_1_GSE138292 \
# -species Cricetulus_griseus \
# -rf /home/ewalsh/FastqProcess/Cricetulus_griseus/GCF_000223135.1_CriGri_1.0_genomic.fna.gz \
# -rg /home/ewalsh/FastqProcess/Cricetulus_griseus/GCF_000223135.1_CriGri_1.0_genomic.gtf.gz \
# -star_args '{"--limitGenomeGenerateRAM": 90000000000}' 

# python3 /home/ewalsh/FastqProcess/main.py STAR -s CHO_1_GSE138292 

# python3 /home/ewalsh/FastqProcess/main.py htseq -s CHO_1_GSE138292 

# python3 /home/ewalsh/FastqProcess/main.py DownstreamAnalysis -s CHO_1_GSE138292 \
# -df BioProject treatment \
# -c treatment insulin control \
# -m /home/ewalsh/FastqProcess/RunFiles/hCMEC_GSE195781/SraRunTable.csv 

# python3 /home/ewalsh/FastqProcess/main.py Normalize -s CHO_1_GSE138292 