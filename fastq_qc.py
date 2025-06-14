#/usr/bin/python3
# Last Update: 
# Capsida Biotherapeutics
# Data Science 
# Evan Walsh (evan.walsh@capsida.com)
# borrowed from /captools/captools-AWS/PyCaptools.py

import os
import sys
import numpy as np
import pandas as pd
import glob
import scipy
from threading import Thread
from matplotlib import pyplot as plt
import pdb
import re
import csv
import threading
import psutil
import time
import shutil
import concurrent
from concurrent.futures import ProcessPoolExecutor

# Save the current working directory
import sys
gen_tools_dir = "/captools/ewalsh/Retrogenix"  # Define the directory where GEN_Tools.py is located
sys.path.insert(0, gen_tools_dir)
import GEN_Tools

def main(args, run_logs):
    """ Run specified program(s). """
    ## command line args 
    indir, root, outroot, study = args.i, args.root, args.o, args.s

    ## initialize Log tools and directories
    my_logs = GEN_Tools.Logs(run_logs)
    breaker = "=" * 120
    my_logs.append_logs(breaker)
    my_logs.append_logs("Initializing fastq QC script")

    ## Initialize inputs and outputs
    if indir is None:
        indir = os.path.join(root, 'expdata', study, "fastq_dump")
    
    fastq_files = GEN_Tools.SearchDir.get_files(parent_dir = indir, extension = "*.fastq.gz", check=None, wd=os.getcwd(), recursive=True)
    
    ## Send fastqc to subdirectory in fastq folders
    outdir = os.path.join(root, 'expdata', study, 'fastqc'); GEN_Tools.SearchDir.check_dir(outdir, critical = False)

    ############################################
    ## Run Fastqc
    path_to_fastqc = shutil.which("fastqc")
    if path_to_fastqc is None:
        sys.exit("Error...fastqc either not installed or not in PATH variable...Exiting")

    threads = []
    for fastq_file in fastq_files:
        cmd = ['fastqc', '-o', outdir, fastq_file]
        t = threading.Thread(target=GEN_Tools.ProcessCall.command_process_string, args=(cmd,))
        t.start()
        threads.append(t)

    # Wait for all threads to complete
    for t in threads:
        t.join()

    ############################################
    ## Run adaptor trimming with fastp

    ## Initialize inputs and outputs
    outdir = os.path.join(root, 'expdata', study, 'fastq_trimmed'); GEN_Tools.SearchDir.check_dir(outdir, critical = False)
    path_to_fastp = shutil.which("fastp")
    if path_to_fastp is None:
        sys.exit("Error...fastp either not installed or not in PATH variable...Exiting")

    subdirs = [name for name in os.listdir(indir) if os.path.isdir(os.path.join(indir, name))]
    srr_ids = [x for x in subdirs if x not in ('fastq_trimmed', 'fastqc')]
    for srr_id in srr_ids:
        ## get SRR ID fastqs
        tmp_dir = os.path.join(outdir, srr_id); GEN_Tools.SearchDir.check_dir(tmp_dir, critical = False)
        fastq_files = GEN_Tools.SearchDir.get_files(parent_dir = os.path.join(indir, srr_id), extension = "*.fastq.gz", check=None, wd=os.getcwd(), recursive=True)
        if len(fastq_files) > 2 and isinstance(fastq_files, list):
            sys.exit("Error... more than 2 fastq files associates with this SRR ID...exiting")
        elif len(fastq_files) == 2:
            fastq_files = sorted(fastq_files) ## sort
            fastq_1 = fastq_files[0]; fastq_trimmed_1 = os.path.join(tmp_dir, os.path.basename(fastq_1))
            fastq_2 = fastq_files[1]; fastq_trimmed_2 = os.path.join(tmp_dir, os.path.basename(fastq_2))
            cmd = ['fastp', "-i", fastq_1, "-I", fastq_2, "-o", fastq_trimmed_1, "-O", fastq_trimmed_2]
            GEN_Tools.ProcessCall.command_process_string(cmd)
        else:
            fastq_file = fastq_files; fastq_trimmed = os.path.join(tmp_dir, os.path.basename(fastq_file))
            cmd = ['fastp', "-i", fastq_file, "-o", fastq_trimmed]
            GEN_Tools.ProcessCall.command_process_string(cmd)

if __name__ == "__main__":
    main()

