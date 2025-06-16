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
import subprocess
import csv
import OS_Tools

class refController:
    """ Class that handles all downloading and creating of reference files onto EC2 VM. """
    def __init__(self, species, outdir, study, ref_fasta, ref_gtf):
        self.species = species
        self.out = outdir 
        self.study = study
        self.ref_fasta = ref_fasta
        self.ref_gtf = ref_gtf

    def get_ref_files(self, data_loc='/home/ewalsh/scratch', ftp_root="https://ftp.ensembl.org/pub/"):
        ## get the path to the species reference files from NCBI ftp site
        targetDir = os.path.join(self.out, self.species)
        OS_Tools.clear_directory(targetDir)

        if self.ref_fasta is None:
            ftp_path = os.path.join(ftp_root, 'current_fasta', self.species.lower(), 'dna/')            
            command = ["wget -r -nd -q -np -A *.dna.primary_assembly.fa.gz -R index.html*,*_from_genomic.fna.gz -e robots=off -P", os.path.join(targetDir, 'fasta'), ftp_path]
            subprocess.run(command)
        else:
            fasta_dir = os.path.join(targetDir, 'fasta')
            OS_Tools.ensure_directory(fasta_dir, critical=False)
            command = ['cp', self.ref_fasta, fasta_dir]
            subprocess.run(command)

        if self.ref_gtf is None:
            # ftp_path = os.path.join(ftp_root, 'genomes', 'refseq', 'vertebrate_mammalian', self.species, 'latest_assembly_versions/')
            ftp_path = os.path.join(ftp_root, 'current_gtf', self.species.lower() + '/')
            command = ["wget -r -nd -q -np -A *.gtf.gz -R index.html*,*chr.gtf.gz,*abinitio.gtf.gz,*scaff.gtf.gz  -e robots=off -P", os.path.join(targetDir, 'gtf'), ftp_path]
            subprocess.run(command)
        else:
            gtf_dir = os.path.join(targetDir, 'gtf')
            OS_Tools.ensure_directory(gtf_dir, critical=False)
            command = ['cp', self.ref_gtf, gtf_dir]
            subprocess.run(command)

    def run_STAR_mkref(self, output_loc, fastaFile, gtfFile, star_args):
        ## delete the genome_dir if it exists
        genome_dir = os.path.join(output_loc, "STAR_index")
        star_tmp_dir = os.path.join(output_loc, '_STARtmp')
        if os.path.exists(genome_dir):
            OS_Tools.clear_directory(genome_dir)

        if os.path.exists(star_tmp_dir):
            OS_Tools.clear_directory(star_tmp_dir)

        OS_Tools.ensure_directory(genome_dir, critical = False)

        cmd = ["STAR", "--runMode", "genomeGenerate", "--genomeDir", genome_dir,
               "--genomeFastaFiles", fastaFile, 
               "--sjdbGTFfile", gtfFile,
               "--sjdbOverhang", str(50),
               "--runThreadN", str(12)]
        

        ## Additional star parameters
        try:
            for key, value in star_args.items():
                cmd.extend([str(key), str(value)])
        except AttributeError:
            pass
        subprocess.run(cmd)

def main(args, run_logs):
    ## Get arguments from argparse
    outroot, study, species, ref_fasta, ref_gtf = args.o, args.s, args.species, args.rf, args.rg
    
    ## initialize Log tools and directories
    my_logs = OS_Tools.Logs(run_logs)
    breaker = "=" * 120

    # outdir = os.path.join(outroot, study)
    OS_Tools.ensure_directory(outroot, critical = False)

    ## create directory to store the reference genome files 
    my_controller = refController(species, outroot, study, ref_fasta, ref_gtf) ## initialize fastq controller
    my_controller.get_ref_files() ## download fastq files from NCBI

    ## Get files
    try:
        searchDir = os.path.join(outroot, args.species, 'fasta', '*.fa.gz');
        refFasta = glob.glob(searchDir)[0]
    except IndexError:
        searchDir = os.path.join(outroot, args.species, 'fasta', '*.fna.gz');
    refFasta = glob.glob(searchDir)[0]
    cmd = ['gunzip', refFasta]; subprocess.run(cmd)
    refFasta = re.sub('.gz','',refFasta)
    searchDir = os.path.join(outroot, args.species, 'gtf', '*.gtf.gz'); 
    refGTF = glob.glob(searchDir)[0]
    cmd = ['gunzip', refGTF]; subprocess.run(cmd)
    refGTF = re.sub('.gz','',refGTF)
    
    #########################################
    ## run spipe in mkref mode
    ## run reference genome generation and send to s3
    ## load the Whole Transcriptome adata object with Scanpy    
    refBuild = Thread(target=my_controller.run_STAR_mkref, args=[outroot, refFasta, refGTF, args.star_args], daemon=True) # daemon means that all threads will exit when the main thread exits
    refBuild.start()
    refBuild.join()

if __name__ == "__main__":
    main()

###