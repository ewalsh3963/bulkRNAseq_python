#/usr/bin/python3
# Last Update: 
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
import warnings
import re
import csv
import subprocess
import gzip
from gtfparse import read_gtf
import io
from pybiomart import Server
import polars as pl
import OS_Tools


class RNASeqNormalizer:
    def __init__(self, counts, gtf_path=None, species = None):
        """
        Args:
            count_matrix_path (str): Path to the count matrix (genes x samples).
            gtf_path (str or None): Path to the GTF file, or None to use external gene lengths.
        """
        self.counts = counts
        self.gtf_path = gtf_path
        self.species = species

        if gtf_path:
            if gtf_path.endswith('.gz'):
                subprocess.run(["gunzip", gtf_path])
                gtf_path = re.sub('.gz$','',gtf_path)
            self.gtf = read_gtf(gtf_path).to_pandas()
            self.gene_lengths = self._get_gene_lengths_from_gtf()
        elif gtf_path is None and species is not None:
            self.gene_lengths = self._get_gene_lengths_external()
        else:
            warnings.warn("Cannot calculte gene lengths. TPM and RPKM nomralizations will fail.")

    def _get_gene_lengths_from_gtf(self):
        """Get gene lengths by summing exon lengths from GTF."""
        exons = self.gtf[self.gtf["feature"] == "exon"]
        vars = [x for x in ['gene_id', 'transcript_id','gene_name'] if x in exons.columns]
        exons = exons[vars + ['start', 'end']]
        exons['length'] = exons['end'] - exons['start'] + 1
        
        gene_lengths = exons.groupby(vars)['length'].sum().reset_index()
        gene_lengths = gene_lengths.groupby(['gene_id'])['length'].median()
        return gene_lengths

    def _get_gene_lengths_external(self):
        """
        Fallback to external gene lengths 
        Retrieves gene lengths using pybiomart for a given species in 'Genus_species' format.

        Args:
            species (str): Species name in 'Genus_species' format (e.g., 'Homo_sapiens').

        Returns:
            pd.DataFrame: DataFrame with gene symbol, Ensembl gene ID, and gene length.
        """

        # Convert to BioMart dataset name format
        species_lower = self.species.lower()
        dataset_name = f"{species_lower.replace('_', '')}_gene_ensembl"

        server = Server(host='http://www.ensembl.org', use_cache=False) #server.list_marts()
        mart = server['ENSEMBL_MART_ENSEMBL'] # mart.list_datasets()
        # Connect to Ensembl BioMart server
        try:
            dataset = mart[dataset_name]
        except KeyError:
            raise ValueError(f"Dataset for species '{self.species}' not found. Please check spelling or supported species.")

        # Query gene info including gene length
        df = dataset.query(attributes=[
            'ensembl_gene_id',
            'external_gene_name',
            'start_position',
            'end_position'
        ])

        # Calculate gene length
        df['gene_length'] = df['Gene end (bp)'] - df['Gene start (bp)'] + 1
        df = df.rename(columns={
            'Gene stable ID': 'ensembl_gene_id',
            'Gene name': 'gene_name'
        })

        return df[['ensembl_gene_id', 'gene_name', 'gene_length']].drop_duplicates()

    def _match_gene_lengths(self):
        """Align gene lengths with count matrix index."""
        return self.gene_lengths.reindex(self.counts.index).fillna(1)

    def to_rpkm(self):
        gene_lengths_kb = self._match_gene_lengths() / 1000
        total_counts_million = self.counts.sum() / 1e6
        rpkm = self.counts.div(gene_lengths_kb, axis=0).div(total_counts_million, axis=1)
        return rpkm

    def to_tpm(self):
        gene_lengths_kb = self._match_gene_lengths() / 1000
        rpk = self.counts.div(gene_lengths_kb, axis=0)
        scaling_factor = rpk.sum()
        tpm = rpk.div(scaling_factor, axis=1) * 1e6
        return tpm

    def to_cpm(self):
        counts_sum = self.counts.sum()
        cpm = self.counts.div(counts_sum, axis=1) * 1e6
        return cpm

def main(args, run_logs):
    ## Get arguments from argparse
    root, outroot, study, indir, gtf_file, format = args.root, args.o, args.s, args.i, args.gtf, args.f

    ## initialize Log tools and directories
    my_logs = OS_Tools.Logs(run_logs)
    breaker = "=" * 120

    if indir is None:
        indir = os.path.join(root, 'PublicRNAseqData', study)

    OS_Tools.ensure_directory(outroot, critical = False)
    outdir = os.path.join(outroot, 'Normalization'); OS_Tools.ensure_directory(outdir, critical = False)

    #########################################
    ## Load, Inspect, and Prepare Data

    ## Join the htseq outputs
    htseq_files = OS_Tools.find_files(parent_dir = indir, 
                                        extension = "htseq_out.txt", 
                                        check=None, 
                                        wd=os.getcwd(), 
                                        recursive=True)
    count_matrix = pd.DataFrame()
    for htseq_file in htseq_files:
        htseq_file.split('/')
        file_parts = htseq_file.split('/')
        srr_id = [item for item in file_parts if item.startswith("SRR")][0]

        dat = pd.read_csv(htseq_file, sep='\t', header = None)
        dat.columns = ['feature', 'count']
        dat['SRR_ID'] = srr_id
        count_matrix = pd.concat([count_matrix, dat])

    count_matrix = count_matrix.pivot_table(values='count', index=['feature'], columns=['SRR_ID']).fillna(0).reset_index()
    count_matrix = count_matrix[~count_matrix['feature'].str.startswith('__')] ## Remove undetermined counts
    count_matrix = count_matrix.set_index('feature')

    if gtf_file is None:
        refDir = os.path.join(root, 'genomes', study); OS_Tools.ensure_directory(refDir, critical = False)
        try:
            gtf_file = OS_Tools.find_files(parent_dir = refDir, extension = "*.gtf.gz", check=None, wd=os.getcwd(), recursive=False)
        except FileNotFoundError:
            try:
                gtf_file = OS_Tools.find_files(parent_dir = refDir, extension = "*.gtf", check=None, wd=os.getcwd(), recursive=False)
            except:
                sys.exit("No GTF found or specified only calcuting CPM normalization...")
                format = "cpm"
    
    normalizer = RNASeqNormalizer(count_matrix, gtf_file)
    if format is None:
        ## Calculate
        cpm_counts = normalizer.to_cpm();
        tpm_counts = normalizer.to_tpm();
        rpkm_counts = normalizer.to_rpkm();

        ## Save 
        outfile = os.path.join(outdir, 'cpm_counts.tsv')
        cpm_counts.to_csv(outfile, sep='\t', header = True, index = True)

        outfile = os.path.join(outdir, 'tpm_counts.tsv')
        tpm_counts.to_csv(outfile, sep='\t', header = True, index = True)

        outfile = os.path.join(outdir, 'rpkm_counts.tsv')
        rpkm_counts.to_csv(outfile, sep='\t', header = True, index = True)

    elif format == "CPM":
        cpm_counts = normalizer.to_cpm();
        outfile = os.path.join(outdir, 'cpm_counts.tsv')
        cpm_counts.to_csv(outfile, sep='\t', header = True, index = True)
    elif format == "TPM":
        tpm_counts = normalizer.to_tpm();
        outfile = os.path.join(outdir, 'tpm_counts.tsv')
        tpm_counts.to_csv(outfile, sep='\t', header = True, index = True)
    elif format == "RPKM":
        rpkm_counts = normalizer.to_rpkm();
        outfile = os.path.join(outdir, 'rpkm_counts.tsv')
        rpkm_counts.to_csv(outfile, sep='\t', header = True, index = True)

if __name__ == "__main__":
    main()

