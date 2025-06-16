#/usr/bin/python3
# Last Update: 
# Capsida Biotherapeutics
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
import shutil
import subprocess
import re
import csv
from threading import Thread, Semaphore
import OS_Tools

def runSTAR(fq, output_dir, genome_dir, star_args):
    
    prefix = re.sub('.fastq.gz','',os.path.basename(fq))
    cmd = ["STAR", 
           "--genomeDir", genome_dir,
           "--runThreadN", "16",
           "--readFilesIn", fq,
           "--outFileNamePrefix", prefix,
           "--outSAMtype", "BAM SortedByCoordinate",
           "--readFilesCommand", "zcat"]
    
    ## Additional star parameters
    try:
        for key, value in star_args.items():
            cmd.extend([key, value])
    except AttributeError:
        pass
        
    subprocess.run(cmd, shell=True, wd=output_dir)

def runSTAR_PE(fqs, output_dir, genome_dir, star_args):
    fqs = sorted(fqs)
    prefix = re.sub('.fastq.gz','',os.path.basename(fqs[0]))
    cmd = ["STAR", 
           "--genomeDir", genome_dir,
           "--runThreadN", "16",
           "--readFilesIn", fqs[0], fqs[1],
           "--outFileNamePrefix", prefix,
           "--outSAMtype", "BAM SortedByCoordinate",
           "--readFilesCommand", "zcat"]
    
    ## Additional star parameters
    try:
        for key, value in star_args.items():
            cmd.extend([key, value])
    except AttributeError:
        pass
        
    subprocess.run(cmd, shell=True, wd=output_dir)

def thread_wrapper(thread_limit, func, args):
    with thread_limit:
        func(*args)  # run the actual target function

def main(args, run_logs):
    ## Get arguments from argparse
    root, outroot, study, indir  = args.root, args.o, args.s, args.i
    
    ## initialize Log tools and directories
    my_logs = OS_Tools.Logs(run_logs)
    breaker = "=" * 120
    my_logs.append_logs("Running STAR alignment...")

    if args.r is None:
        refDir = os.path.join(root, 'genomes', study)
        my_logs.append_logs("Reference directory for alignment:\n{}".format(refDir))
        genome_dir = os.path.join(refDir, 'STAR_index') ## if not get it from the default location
    else: 
        if args.r.startswith('s3'): ## pull from S3 if needed
            refDir = os.path.join(root, 'genomes', study)
            my_logs.append_logs("Reference directory for alignment:\n{}".format(refDir))
            cmd = ['aws', 's3', 'sync', args.r, refDir, '--no-progress']
            t = Thread(target=subprocess.run, args=[cmd], daemon=True) # daemon means that all threads will exit when the main thread exits
            t.start()
            t.join()

            genome_dir = os.path.join(refDir,'STAR_index')
        else: 
            refDir = args.r
            my_logs.append_logs("Reference directory for alignment:\n{}".format(refDir))
            genome_dir = os.path.join(refDir,'STAR_index')

    #########################################
    ## Check to make sure STAR and samtools are installed and in PATH 
    path_to_STAR = shutil.which("STAR")
    if path_to_STAR is None:
        sys.exit("Error...STAR either not installed or not in PATH variable...Exiting")

    path_to_samtools = shutil.which("samtools")
    if path_to_samtools is None:
        sys.exit("Error...STAR either not installed or not in PATH variable...Exiting")

    #########################################
    ## Initialize inputs and outputs
    OS_Tools.ensure_directory(outroot, critical = False)
    if indir is None:
        indir = os.path.join(root, 'expdata', study, "fastq_trimmed")
    
    #########################################
    ## Run alignment
    thread_limit = Semaphore(2)
    threads = []
    srr_ids = [name for name in os.listdir(indir) if os.path.isdir(os.path.join(indir, name))]
    for srr_id in srr_ids:
        fqs = OS_Tools.find_files(
            parent_dir=os.path.join(indir, srr_id),
            extension="*.fastq.gz",
            check=None,
            wd=os.getcwd(),
            recursive=True
        )
        outdir = os.path.join(outroot, srr_id)
        OS_Tools.ensure_directory(outdir, critical=False)
        if len(fqs) > 1 and isinstance(fqs, list):
            fqs = sorted(fqs)
            t = Thread(target=thread_wrapper, args=(thread_limit, runSTAR_PE, [fqs, outdir, genome_dir, args.star_args]), daemon=True)
        else:
            t = Thread(target=thread_wrapper, args=(thread_limit, runSTAR, [fqs, outdir, genome_dir, args.star_args]), daemon=True)

        t.start()
        threads.append(t)

    ## Wait for all threads to complete
    for t in threads:
        t.join()

    ## Index Bam Files after alignment
    for srr_id in srr_ids:
        outdir = os.path.join(outroot, srr_id);
        bam_file = OS_Tools.find_files(parent_dir = outdir, extension = "*.bam", check=None, wd=os.getcwd())
        cmd = ["samtools", "index", bam_file]
        subprocess.run(cmd, shell=True, wd=outdir)



if __name__ == "__main__":
    main()

#   nohup /usr/bin/STAR --genomeDir /home/ewalsh/scratch/genomes/${study}_TargetGene/STAR_index \
#   --runThreadN 16 --readFilesIn ${fq} \
#   --outFileNamePrefix ${sublib} \
#   --outSAMtype BAM SortedByCoordinate \
#   --readFilesCommand zcat \
#   --outFilterMultimapNmax 3 \
#   --outSAMunmapped Within &