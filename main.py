#!/usr/bin/env python3
import sys
import argparse as ap
import pandas as pd
import numpy as np
import subprocess
import plot
import motif_enrichment as enrichment
from sklearn.preprocessing import StandardScaler


def parse_args():
    '''
    Argument Parser
    '''
    parser = ap.ArgumentParser(description="correct way to parse",
                               prog='time_series')

    parser.add_argument('-o',
                        '--out_dir',
                        type=str,
                        help="Output directory",
                        default='./out',
                        required=False)

    parser.add_argument('-c',
                        '--counts',
                        type=str,
                        help="Input read counts filename",
                        default='./data/raw_counts.txt',
                        required=False)

    parser.add_argument('-k',
                        '--kmeans',
                        type=str,
                        help="Optional: Manually select a k for clustering",
                        default=None,
                        required=False)

    parser.add_argument('-n',
                        '--num_genes',
                        type=str,
                        help="Number of genes to plot initial heatmap etc.",
                        default=500,
                        required=False)

    parser.add_argument('-a',
                        '--ame',
                        type=str2bool,
                        help="Enable AME after installing AME <True/False>",
                        default=False,
                        required=False)

    return parser.parse_args()


def str2bool(v):
    '''
    small function for boolean input
    '''
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def read_counts(countfile):
    '''
    Read in the raw read counts from featureCounts
    '''
    counts_header = None
    counts = []
    for l in open(countfile):
        if counts_header is None:
            counts_header = l.rstrip().split("\t")
        else:
            gene_count = l.rstrip().split("\t")
            gene_count[2] = gene_count[2].split(', ', 1)[0]
            # Remove version number from refseq accession number
            counts.append(gene_count)
    return counts_header, counts


def sort_counts(in_counts):
    '''
    sort counts by variance across time
    '''
    # Turning off pandas SettingWithCopyWarning
    pd.set_option('mode.chained_assignment', None)
    # Some array's SNR is too low, filter and resort
    order = ['ACCNUM', '0.15H.IFN_1', '0.5H.IFN_1', '0.75H.IFN_1',
             '1H.IFN_1', '1.25H.IFN_1', '1.5H.IFN_1', '2.25H.IFN_1',
             '2.75H.IFN_1', '3H.IFN_1', '3.5H.IFN_1', '5H.IFN_1',
             '5.5H.IFN_1', '6H.IFN_1', '6.5H.IFN_1', '7H.IFN_1',
             '8H.IFN_1', '9H.IFN_1', '10H.IFN_1', '11H.IFN_1',
             '12H.IFN_1', '13H.IFN_1', '14H.IFN_1', '15H.IFN_1']
    in_counts['SYMBOL'] = in_counts['SYMBOL'].replace([' ', ',', ';'],
                                                      '_', regex=True)
    in_counts = in_counts[order]
    # Sort genes by variance
    x = in_counts.iloc[:, 1:].values
    x = StandardScaler().fit_transform(x)
    x_sub = np.subtract(x, x[:, 0].reshape((len(x), 1)))
    x_var = np.var(x_sub, 1, dtype=np.float64)
    in_counts['Variance'] = x_var
    in_counts.sort_values(by=['Variance'], inplace=True, ascending=False)
    return in_counts.iloc[0:2000, :-2]


def main():
    args = parse_args()
    header, counts = read_counts(args.counts)
    counts = pd.DataFrame(counts, columns=header)
    counts = sort_counts(counts)

    num_genes = int(args.num_genes)-1
    out_dir = args.out_dir

    # plot heatmap
    print('Generating Gene Expression Heatmap ...')
    plot.plot_heatmap(out_dir, counts, num_genes)

    # Plot gene expression tragectory
    print('Generating Gene Expression Trajectory ...')
    plot.plot_trajectory(out_dir, counts, num_genes)

    # Save formatted and sorted count file
    counts.to_csv('data/counts_clust.txt', index=False, sep='\t')

    # Call clust to cluster the genes
    print('Clustering Genes (this might take awhile) ...')
    subprocess.call(['mkdir', 'clust_out'])
    if args.kmeans is not None:
        subprocess.call(['clust', 'data/counts_clust.txt',
                         '-o', './clust_out', '-K', args.kmeans])
    else:
        subprocess.call(['clust', 'data/counts_clust.txt',
                        '-o', './clust_out'])

    # Make sure everythin is all set for motif enrichment
    clust_out = './clust_out/Clusters_Objects.tsv'
    ref_bed = './ref/mm10.refseq.bed'
    ref_fa = './ref/mm10.fa'
    ref_motif = './ref/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme'

    try:
        open(clust_out)
    except FileNotFoundError:
        print('Clustering failed, '
              + 'Clusters_Objects.tsv not found in the output directory')

    print('Motif Enrichment Analysis (this might take awhile) ...')
    enrichment.run_motif_enrichment(clust_out,
                                    ref_bed,
                                    ref_fa, ref_motif,
                                    args.ame,
                                    100, 100)


if __name__ == '__main__':
    main()
