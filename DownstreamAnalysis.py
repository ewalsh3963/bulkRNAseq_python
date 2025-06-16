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
from mygene import MyGeneInfo
from threading import Thread
from matplotlib import pyplot as plt
import scanpy as sc
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import seaborn as sns
from kneed import KneeLocator
from pydeseq2.dds import DeseqDataSet  # correct for version >= 0.4.0
from pydeseq2.ds import DeseqStats
import gseapy as gp
from gseapy.plot import gseaplot
import sys



def main(args, run_logs):
    ## Get arguments from argparse
    root, outroot, study, indir, min_samples, metadata_file  = args.root, args.o, args.s, args.i, args.ms, args.m
    
    #########################################
    ## Initialize inputs and outputs
    if indir is None:
        indir = os.path.join(root, 'PublicRNAseqData', study)

    OS_Tools.ensure_directory(outroot, critical = False)
    outdir = os.path.join(outroot, 'DownstreamAnalysis'); OS_Tools.ensure_directory(outdir, critical = False)

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

    #########################################
    ## Exploratory Data Analysis

    # Drop the Gene Name column for counting
    countlist_no_name = count_matrix.iloc[:, 1:]

    # Calculate total counts and log transform
    total_counts = countlist_no_name.sum(axis=0)
    log_counts = countlist_no_name.apply(lambda x: np.log2(x + 1))

    # Create main visualization figure
    fig1, axes = plt.subplots(2, 2, figsize=[8,8])

    # Panel 1: Total counts per sample
    sns.barplot(x=countlist_no_name.columns, y=total_counts, color='skyblue', ax=axes[0,0])
    axes[0,0].set_ylabel('Total Counts')
    axes[0,0].set_title('Total Counts per Sample')
    axes[0,0].tick_params(axis='x', rotation=85)

    # Panel 2: Log transformed counts distribution
    log_counts.boxplot(ax=axes[0,1])
    axes[0,1].set_ylabel('Log2(Counts + 1)')
    axes[0,1].set_title('Log Transformed Counts per Sample')
    axes[0,1].tick_params(axis='x', rotation=85)

    # Panel 3: Sample correlation heatmap
    correlation_matrix = log_counts.corr()
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0.5, vmin=0, vmax=1, ax=axes[1,0])
    axes[1,0].set_yticklabels(axes[1, 0].get_yticklabels(), rotation=45)
    axes[1,0].set_title('Sample Correlation Matrix')

    # Panel 4: PCA plot
    pca = PCA(n_components=2)
    scaler = StandardScaler()
    pca_result = pca.fit_transform(scaler.fit_transform(log_counts.T))
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=log_counts.columns)
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', s=100, ax=axes[1,1])
    for idx, row in pca_df.iterrows():
        axes[1,1].annotate(idx, (row['PC1'], row['PC2']))
    axes[1,1].set_title(f'PCA Plot\nPC1 ({pca.explained_variance_ratio_[0]:.1%}) vs PC2 ({pca.explained_variance_ratio_[1]:.1%})')
    plt.tight_layout()

    outfile = os.path.join(outdir, 'SampleStats.pdf')
    plt.savefig(outfile)

    # Create dendrogram figure
    fig2 = plt.figure(figsize=(8, 6))
    h_clustering = linkage(log_counts.T, 'ward')
    dendrogram(h_clustering, labels=countlist_no_name.columns)
    plt.xticks(rotation=90)
    plt.ylabel('Distance')
    plt.title('Sample Clustering Dendrogram')
    outfile = os.path.join(outdir, 'SampleHierarchy.pdf')
    plt.savefig(outfile)
    
    # Generate QC metrics and save 
    qc_stats = {
        'total_reads': total_counts.sum(),
        'mean_reads_per_sample': total_counts.mean(),
        'cv_reads': total_counts.std() / total_counts.mean(),
        'min_sample_correlation': correlation_matrix.min().min(),
        'max_sample_correlation': correlation_matrix.max().min(),
        'pc1_variance': pca.explained_variance_ratio_[0],
        'pc2_variance': pca.explained_variance_ratio_[1]}
    
    qcDat =  pd.DataFrame.from_dict(qc_stats, orient='index').reset_index()
    qcDat.columns = ['Stat','Value']
    outfile = os.path.join(outdir, 'qcStats.csv')
    qcDat.to_csv(outfile, sep=',', header = True)

    #########################################
    ## Quality Control, Filtering, and Normalization

    # convert raw counts to CPM to normalize the data
    cpm = countlist_no_name.apply(lambda x: (x / x.sum()) * 1e6) #convert raw counts to CPM to normalize
    # define a range of CPM thresholds to test, from 0 to 5 with increments of 0.1
    thresholds = np.arange(0, 5, 0.1)
    # initialize list to store the # of genes retained for ea/ threshold
    genes_retained = []

    # loop through ea/ threshold value to determine the # of genes retained
    for min_cpm in thresholds:
        # create mask where CPM > min_cpm in at least min_samples samples
        mask = (cpm > min_cpm).sum(axis=1) >= min_samples
        # count # of genes that meet the criteria and append to the list
        genes_retained.append(mask.sum())

    # Use KneeLocator to find the inflection point
    kneedle = KneeLocator(thresholds, genes_retained, curve='convex', direction='decreasing')
    knee_cpm = kneedle.knee
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, genes_retained, marker='o', color='green', label='Genes retained')
    plt.axvline(x=1.0, color='red', linestyle='--', label='CPM = 1')
    if knee_cpm is not None:
        plt.axvline(x=knee_cpm, color='blue', linestyle='--', label=f'Knee CPM = {knee_cpm:.2f}')
        plt.scatter(knee_cpm, kneedle.knee_y, color='blue', zorder=5)

    plt.xlabel('Threshold (CPM)')
    plt.ylabel('Num Genes Retained')
    plt.title('Genes Retained vs. CPM Threshold')
    plt.legend()
    plt.tight_layout()
    outfile = os.path.join(outdir, 'CPM_Threshold.pdf')
    plt.savefig(outfile)

    #########################################
    ## DESeq2 Analysis
    filtered_counts = count_matrix[mask].set_index('feature').T
    if metadata_file is None:
        return
    else:
        metadata = pd.read_csv(metadata_file, sep=',', header = 0)
        metadata = metadata.set_index('Run')

        ## set design factors as categorical and remove 
        design_factors = args.df
        for df in design_factors:
            metadata[df] = metadata[df].astype('category')
            if len(metadata[df].unique()) < 2:
                design_factors.remove(df)
        
        ## Create DDS object
        dds = DeseqDataSet(counts=filtered_counts, metadata = metadata, design_factors = design_factors)

        ## Run DESeq2
        dds.deseq2()
        contrast = tuple(args.contrast)
        stat_res = DeseqStats(dds, contrast=contrast)

        ## save results 
        # summary_dat = stat_res.results_df  
        # outfile = os.path.join(outdir, 'DEseq2_summary.csv')
        # summary_dat.to_csv(outfile, sep=',', header = True)
        
        stat_res.summary()
        res = stat_res.results_df
        outfile = os.path.join(outdir, 'DEseq2_results.csv')
        res.to_csv(outfile, sep=',', header = True)

        ## Significantly differentially expressed genes
        padj_threshold, abs_logFC_threshold = args.padj_threshold, args.abs_logFC
        res['abs_log2fc'] = res['log2FoldChange'].abs()
        res['significant'] = np.where(
            (res['padj'] <= padj_threshold) & (res['abs_log2fc'] > abs_logFC_threshold),
            True,
            False
        )
        sigs = res[res['significant']]
        outfile = os.path.join(outdir, 'DEseq2_sig_results.csv')
        sigs.to_csv(outfile, sep=',', header = True)

        #########################################
        ## Transform count data 
        sc.tl.pca(dds)
        sc.pl.pca(dds, color = contrast[0], size = 200)
        outfile = os.path.join(outdir, 'DDS_PCA.pdf')
        plt.tight_layout()
        plt.savefig(outfile)
        plt.close()

        #########################################
        ## Gene set enrichment analysis (GSEA)
        ranking = res[['stat']].dropna().sort_values('stat', ascending = False)

        # your index should be the Ensembl gene IDs
        mg = MyGeneInfo()
        ensembl_ids = ranking.index.tolist()

        # Query Ensembl to Symbol
        query_result = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')

        # Build a mapping dict from results
        id_to_symbol = {item['query']: item.get('symbol', None) for item in query_result if 'symbol' in item}

        # Filter and rename
        ranking = ranking[ranking.index.isin(id_to_symbol.keys())]
        ranking.index = ranking.index.map(id_to_symbol)
        ranking = ranking[~ranking.index.duplicated()]  # remove duplicates
        enr_res = gp.prerank(rnk = ranking, gene_sets = ['GO_Biological_Process_2021'], seed = 6, permutation_num = 100)
        
        out = []
        for term in list(enr_res.results):
            out.append([term,
                    enr_res.results[term]['fdr'],
                    enr_res.results[term]['es'],
                    enr_res.results[term]['nes']])

        out_df = pd.DataFrame(out, columns = ['Term','fdr', 'es', 'nes']).sort_values('fdr').reset_index(drop = True)
        outfile = os.path.join(outdir, 'GSEA.csv')
        out_df.to_csv(outfile, sep=',', header = True)

        #########################################
        ## Visualize Differential Expression
        
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        scatter_params = {'alpha': 0.8,'edgecolor': None,'palette': 'viridis'}

        # Panel 1: Global Expression Landscape (Volcano Plot)
        sns.scatterplot(data=res, x='log2FoldChange', y='padj', hue='log2FoldChange',ax=axes[0,0], **scatter_params)
        axes[0,0].axhline(y=padj_threshold, color='red', linestyle='--', linewidth=1)
        axes[0,0].axvline(x=abs_logFC_threshold, color='blue', linestyle='--', linewidth=1)
        axes[0,0].axvline(x=-abs_logFC_threshold, color='blue', linestyle='--', linewidth=1)
        axes[0,0].set_xlabel('log2 Fold Change')
        axes[0,0].set_ylabel('Adjusted P-value')
        axes[0,0].set_title('Global Expression Landscape')

        # Panel 2: Fold Change Distribution (All Genes)
        sns.histplot(data=res, x='abs_log2fc',bins=50, kde=True,ax=axes[0,1])

        # Add vertical line at fold change threshold
        axes[0,1].axvline(x=abs_logFC_threshold, color='red', linestyle='--', linewidth=1)

        axes[0,1].set_title('Distribution of Absolute log2FC (All Genes)')
        axes[0,1].set_xlabel('Absolute log2 Fold Change')
        axes[0,1].set_ylabel('Gene Frequency')

        # Panel 3: MA Plot
        # pdb.set_trace()
        # res['mean_expression'] = np.log2((res['mean_treated'] + res['mean_control'])/2 + 1)
        sns.scatterplot(data=res, x='baseMean', y='log2FoldChange', hue='significant' if 'significant' in res.columns else None, ax=axes[1,0], **scatter_params)
        axes[1,0].axhline(y=0, color='red', linestyle='--', linewidth=1)
        axes[1,0].set_title('MA Plot (Mean vs Fold Change)')
        axes[1,0].set_xlabel('Mean Expression (log2)')
        axes[1,0].set_ylabel('log2 Fold Change')

        # Panel 4: Distribution of Adjusted P-values
        sns.histplot(data=res ,x='padj',bins=50, kde=True, ax=axes[1,1])

        # Add vertical line at significance threshold
        axes[1,1].axvline(x=padj_threshold, color='red', linestyle='--', linewidth=1)
        axes[1,1].set_title('Distribution of Adjusted P-values')
        axes[1,1].set_xlabel('Adjusted P-value')
        axes[1,1].set_ylabel('Gene Frequency')

        plt.tight_layout()
        outfile = os.path.join(outdir, 'DESeq2_results.pdf')
        plt.tight_layout()
        plt.savefig(outfile)
        plt.close()

        # Generate comprehensive analytical metrics
        summary_stats = {
            'total_genes': len(res),
            'significant_genes': len(sigs),
            'mean_fold_change_all': res['abs_log2fc'].mean(),
            'median_fold_change_all': res['abs_log2fc'].median(),
            'max_fold_change': res['abs_log2fc'].max(),
            'mean_fold_change_sig': sigs['abs_log2fc'].mean(),
            'median_padj': res['padj'].median(),
            'genes_below_alpha': sum(res['padj'] < padj_threshold)}


        summary_stats =  pd.DataFrame.from_dict(summary_stats, orient='index').reset_index()
        summary_stats.columns = ['Stat','Value']
        outfile = os.path.join(outdir, 'summaryStats.csv')
        summary_stats.to_csv(outfile, sep=',', header = True)


if __name__ == "__main__":
    main()
