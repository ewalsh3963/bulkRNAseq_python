#/usr/bin/python3
# Data Science 
# Evan Walsh (evanwalsh396@gmail.com)

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
import psutil
import time
import subprocess
import OS_Tools
import concurrent
from concurrent.futures import ProcessPoolExecutor

class Download:
    def __init__(self, srr_table):
        self.srr_table = srr_table
    
    def prefetch(self, srr_id, data_loc, my_logs):
        ## run prefetch to download a sra file before creating fastq
        cmd = ["prefetch", "--max-size 100G", '-O', data_loc, srr_id]
        my_logs.append_logs("Prefetch command:\n")
        my_logs.append_logs(cmd)
        subprocess.run(cmd)
        
        ## move the SRR ID folder 
        outdir = os.path.join(data_loc, srr_id)
        sra_dir = outdir + '_sra'
        cmd = ['mv', outdir, sra_dir]
        subprocess.run(cmd)

        ## run fastq -dump
        sra_file = os.path.join(sra_dir, srr_id + '.sra')
        cmd = ["fastq_dump", "-t", "/home/ewalsh/scratch/tmp", " --include-technical", '-O', outdir, '--split-files', sra_file]
        my_logs.append_logs("Fastq-dump command:\n")
        cmd_str = ' '.join([str(x) for x in cmd])
        my_logs.append_logs(cmd_str)
        subprocess.run(cmd)

        ## Zip fastq files
        cmd = ['gzip', '-f', os.path.join(outdir, '*.fastq')]
        cmd = ' '.join(cmd)
        subprocess.run(cmd, shell = True, cwd=outdir)

    def fastq_dump(self, srr_id, data_loc, my_logs, N=None, X=None):
        ## Make SRR ID Directory
        outdir = os.path.join(data_loc, srr_id); OS_Tools.ensure_directory(outdir, critical = False)
        cmd = ["fastq-dump", "-N", N, "-X", X, '-O', outdir, '--split-files', srr_id]
        my_logs.append_logs("Fastq-dump command:\n")
        cmd_str = ' '.join([str(x) for x in cmd])
        my_logs.append_logs(cmd_str)
        subprocess.run(cmd)

        ## Zip fastq files
        cmd = ['gzip', '-f', os.path.join(outdir, '*.fastq')]
        cmd = ' '.join(cmd)
        subprocess.run(cmd, shell = True, cwd=outdir)

    def sra_toolkit(self, data_loc='~/rnaseq_tools'):
        ## download the sra-toolkit
        outfile = os.path.join(data_loc, "sratoolkit.tar.gz")
        cmd = ["wget", "-O", outfile, "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz"]
        subprocess.run(cmd)

        ## untar the zipped folder 
        cmd = ["tar", "-vxzf",  data_loc + '/' + "sratoolkit.tar.gz", "-C", data_loc]
        subprocess.run(cmd)

        ## Add to path
        toolkit_path = os.path.abspath(data_loc + '/' + "sratoolkit.3.0.7-ubuntu64/bin")

        # Add it to the PATH environment variable
        os.environ["PATH"] = os.environ["PATH"] + os.pathsep + toolkit_path

def main(args, run_logs):
    """ Run specified program(s). """
    ## command line args 
    outroot, study, srr_file, N, X = args.o, args.s, args.srr, args.N, args.X
    outdir = os.path.join(outroot, 'fastq_dump'); OS_Tools.ensure_directory(outdir, critical = False)
    ## initialize Log tools and directories
    my_logs = OS_Tools.Logs(run_logs)
    breaker = "=" * 120
    my_logs.append_logs(breaker)
    my_logs.append_logs("Initializing fastq-dump script")

    Downloader = Download(srr_file)

    ## download SRA-toolkit if not already
    sra_path = "sratoolkit.3.2.1-ubuntu64"
    is_in_path = any([sra_path in x for x in os.environ["PATH"].split(os.pathsep)])
    if not is_in_path:
        Downloader.sra_toolkit()

    ## download fastq files with fastq dump 
    srr_runs = pd.read_csv(srr_file, sep=',', header = 0, index_col = False, usecols=['Run'], comment='#')
    srr_runs = srr_runs[srr_runs['Run'].notna()]
    srr_runs = srr_runs['Run'].tolist()
    if N is None and X is None:
        with ProcessPoolExecutor(5) as executor:
            tasks = [executor.submit(Downloader.prefetch, srr_id, outdir, my_logs) for srr_id in srr_runs]
            for t in tasks:
                t.result()  # wait for all to finish, will also raise exceptions if any
    else:
        with ProcessPoolExecutor(5) as executor:
            tasks = [executor.submit(Downloader.fastq_dump, srr_id, outdir, my_logs, N, X) for srr_id in srr_runs]
            for t in tasks:
                t.result()  # wait for all to finish, will also raise exceptions if any
    
if __name__ == "__main__":
    main()

