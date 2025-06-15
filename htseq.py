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
import subprocess
import re
import shutil
import csv
import OS_Tools




def runHTSEQ(bam_file, gtf_file, outfile, stranded = "no", htseq_args = None):
    if htseq_args is None:
        cmd = ["htseq-count", 
            "-f", "bam", # Specifies the format of the input alignment file; in this case, it's BAM (Binary Alignment/Map format). 
            "-r", "pos", # Specifies the sorting order of the input file; "pos" means the reads are sorted by position.
            "-s", stranded, # Specifies whether the data is strand-specific:
            bam_file,
            gtf_file, ">", outfile]
    else:
        cmd = ["htseq-count"]
        # Add user-specified HTSeq arguments from JSON (e.g., {"-f": "bam", "-r": "pos", "-s": "no"})
        for key, value in htseq_args.items():
            cmd.extend([key, value])
        
        # Add BAM and GTF files
        cmd.extend([bam_file, gtf_file, ">", outfile])

    subprocess.run(cmd, shell=True)

def main(args, run_logs):
    ## Get arguments from argparse
    root, outroot, study, indir, gtf_file, bed_file, htseq_args = args.root, args.o, args.s, args.i, args.gtf, args.bed, args.htseq_args
    
    ## initialize Log tools and directories
    my_logs = OS_Tools.Logs(run_logs)
    breaker = "=" * 120

    #########################################
    ## Check to make sure STAR and samtools are installed and in PATH 
    path_to_STAR = shutil.which("STAR")
    if path_to_STAR is None:
        sys.exit("Error...STAR either not installed or not in PATH variable...Exiting")

    #########################################
    ## Initialize inputs and outputs
    OS_Tools.ensure_directory(outroot, critical = False)
    if indir is None:
        indir = os.path.join(root, 'results', study)

    #########################################
    ## run htseq-count
    threads = []
    srr_ids = [name for name in os.listdir(indir) if os.path.isdir(os.path.join(indir, name))]
    for srr_id in srr_ids:
        bam = OS_Tools.find_files(parent_dir = os.path.join(indir, srr_id), extension = "*.bam", check=None, wd=os.getcwd(), recursive=True)
        outdir = os.path.join(outroot, srr_id); OS_Tools.ensure_directory(outdir, critical = False)
        outfile = os.path.join(outdir, "htseq_out.txt")

        refDir = os.path.join(root, 'genomes', study); OS_Tools.ensure_directory(refDir, critical = False)
        ## need a bed12 file to determine if RNAseq data is stranded or not 
        if gtf_file is None and bed_file is None:
            ## Get GTF file 
            tmp_gtf = "/scratch/tmp.gtf"
            try:
                gtf_file = OS_Tools.find_files(parent_dir = refDir, extension = "*.gtf.gz", check=None, wd=os.getcwd(), recursive=False)

                ## Unzip it and store in scratch
                cmd = ['gunzip',
                    '-c', ## keep original
                    gtf_file, ">", tmp_gtf]
                subprocess.run(cmd, shell=True)
            except FileNotFoundError:
                gtf_file = OS_Tools.find_files(parent_dir = refDir, extension = "*.gtf", check=None, wd=os.getcwd(), recursive=False)
                cmd = ['cp',gtf_file, tmp_gtf]
                subprocess.run(cmd, shell=True)
                
            ## define bed file 
            bed_name = re.sub('.gtf.*','.bed',os.path.basename(gtf_file))
            bed_file = os.path.join(refDir, bed_name)

            ## Convert gtf file to bed12 format with Bedops
            path_to_gtf2bed = shutil.which("gtf2bed")
            if path_to_gtf2bed is None:
                sys.exit("Error...gtf2bed either not installed or not in PATH variable...Exiting")
            
            ## Run gtf2bed
            cmd = ["/usr/local/bin/gtf2bed", "<", tmp_gtf, ">", bed_file]
            subprocess.run(cmd, shell=True)

            ## Infer experiment 
            cmd = ["infer_experiment.py", "-r", bed_file, "-i", bam, "> /scratch/tmp.out"]
            subprocess.run(cmd, shell=True)

            with open("/scratch/tmp.out", "r") as file:
                lines = [line.strip() for line in file if line.strip()]
                stranded_fraction = float(lines[3].split(':')[1])
                if stranded_fraction >= 0.9:
                    stranded = 'yes'
                else: 
                    stranded = 'no'
        elif gtf_file is not None and bed_file is None:
            ## define bed file 
            bed_name = re.sub('.gtf.*','.bed',os.path.basename(gtf_file))
            bed_file = os.path.join(refDir, bed_name)

            ## Convert gtf file to bed12 format with Bedops
            path_to_gtf2bed = shutil.which("gtf2bed")
            if path_to_gtf2bed is None:
                sys.exit("Error...gtf2bed either not installed or not in PATH variable...Exiting")

            ## Run gtf2bed
            cmd = ["/usr/local/bin/gtf2bed", "<", gtf_file, ">", bed_file]
            subprocess.run(cmd, shell=True)

            ## Infer experiment 
            cmd = ["infer_experiment.py", "-r", bed_file, "-i", bam, "> /scratch/tmp.out"]
            subprocess.run(cmd, shell=True)
            with open("/scratch/tmp.out", "r") as file:
                lines = [line.strip() for line in file if line.strip()]
                stranded_fraction = float(lines[3].split(':')[1])
                if stranded_fraction >= 0.9:
                    stranded = 'yes'
                else: 
                    stranded = 'no'
        else:
                ## Infer experiment 
                cmd = ["infer_experiment.py", "-r", bed_file, "-i", bam, "> /scratch/tmp.out"]
                subprocess.run(cmd, shell=True)
                with open("/scratch/tmp.out", "r") as file:
                    lines = [line.strip() for line in file if line.strip()]
                    stranded_fraction = float(lines[3].split(':')[1])
                    if stranded_fraction >= 0.9:
                        stranded = 'yes'
                    else: 
                        stranded = 'no'
        
        if htseq_args == {}:  ## If htseq CL arguments are not specified determine if the experiment is stranded  
            t = Thread(target=runHTSEQ, args=[bam, gtf_file, outfile, stranded, htseq_args], daemon=True) # daemon means that all threads will exit when the main thread exits
            t.start()
            threads.append(t)

        else: ## otherwise run with arguments specifed at CL
            t = Thread(target=runHTSEQ, args=[bam, gtf_file, outfile, "no", htseq_args], daemon=True) # daemon means that all threads will exit when the main thread exits
            t.start()
            threads.append(t)

    # Wait for all threads to complete
    for t in threads:
        t.join()

if __name__ == "__main__":
    main()
