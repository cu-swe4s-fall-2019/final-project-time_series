#!/usr/bin/env python3
import argparse as ap
import pandas as pd
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
sns.set()


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

    parser.add_argument('-n',
                        '--num_genes',
                        type=str,
                        help="Number of genes to plot",
                        default=500,
                        required=False)

    return parser.parse_args()


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
            # Remove version number from refseq accession number
            # gene_count[0] = gene_count[0].split('.', 1)[0]
            counts.append(gene_count)
    return counts_header, counts


def sort_counts(in_counts):
    '''
    sort counts by variance across time
    '''
    # Turning off pandas SettingWithCopyWarning
    pd.set_option('mode.chained_assignment', None)
    # Some array's SNR is too low, filter and resort
    order = ['SYMBOL', '0.15H.IFN_1', '0.5H.IFN_1', '0.75H.IFN_1',
             '1H.IFN_1', '1.25H.IFN_1', '1.5H.IFN_1', '2.25H.IFN_1',
             '2.75H.IFN_1',
             '3H.IFN_1', '3.5H.IFN_1', '5H.IFN_1', '5.5H.IFN_1', '6H.IFN_1',
             '6.5H.IFN_1', '7H.IFN_1', '8H.IFN_1', '9H.IFN_1', '10H.IFN_1',
             '11H.IFN_1', '12H.IFN_1', '13H.IFN_1', '14H.IFN_1', '15H.IFN_1']
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
    return in_counts.iloc[:, :-2]


def plot_pca(out, counts):
    '''
    Function to plot PCA for microarray datasets
    '''
    counts = counts.iloc[:, 1:]
    x = np.transpose(counts.values)
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents,
                               columns=['PC1', 'PC2'])
    finalDf = pd.concat([principalDf, pd.DataFrame(counts.columns)], axis=1)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('PC1', fontsize=15)
    ax.set_ylabel('PC2', fontsize=15)
    ax.set_title('2 Component PCA', fontsize=20)
    ax.scatter(finalDf.loc[:, 'PC1'], finalDf.loc[:, 'PC2'],
               s=50)
    for i, txt in enumerate(counts.columns):
        ax.annotate(txt, (finalDf.loc[i, 'PC1'], finalDf.loc[i, 'PC2']))

    ax.grid()
    plt.savefig(out+'/sample_pca.png', bbox_inches='tight')
    plt.close()


def plot_heatmap(out, counts, num_genes):
    '''
    Plot a heatmap of gene expression over time for top num_genes genes with
    the greatest variance
    '''
    select_counts = counts.iloc[0:num_genes, 1:].astype(float)
    select_counts = select_counts.sub(select_counts.iloc[:, 0], axis='rows')
    gene_names = counts.iloc[0:num_genes, 1]
    # Normalize it by row:
    # (not sure if it is the best way,
    # please feel free to give me a better method.)
    # 1: substract mean
    df_norm_row = select_counts.sub(select_counts.mean(axis=1), axis=0)
    # 2: divide by standard dev
    df_norm_row = df_norm_row.div(select_counts.std(axis=1), axis=0)

    xlabel = [x.replace('.IFN_1', '') for x in select_counts.columns]
    xlabel = [x.replace('H', ' hr') for x in xlabel]
    if num_genes > 50:
        ax = sns.clustermap(df_norm_row,
                            metric="correlation", col_cluster=False,
                            xticklabels=xlabel,
                            cmap="RdBu_r")
    else:
        ax = sns.clustermap(df_norm_row,
                            metric="correlation", col_cluster=False,
                            xticklabels=xlabel,
                            yticklabels=gene_names,
                            cmap="RdBu_r")
    plt.savefig(out+'/gene_heatmap.png', bbox_inches='tight')
    plt.close()


def plot_trajectory(out, counts, num_genes):
    select_counts = counts.iloc[0:num_genes, 1:].astype(float)
    xlabel = [x.replace('.IFN_1', '') for x in select_counts.columns]
    xlabel = [float(x.replace('H', '')) for x in xlabel]
    select_counts = select_counts.sub(select_counts.iloc[:, 0], axis='rows')
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylabel('Gene Expression Fold Change', fontsize=15)
    ax.set_xlabel('Time (hrs)', fontsize=15)
    ax.plot(xlabel, select_counts.T)
    plt.savefig(out+'/gene_trajectory.png', bbox_inches='tight')
    plt.close()


def main():
    args = parse_args()
    header, counts = read_counts(args.counts)
    counts = pd.DataFrame(counts, columns=header)
    counts = sort_counts(counts)

    num_genes = int(args.num_genes)-1
    out_dir = args.out_dir
    # plot heatmap
    print('Generating Gene Expression Heatmap ...')
    plot_heatmap(out_dir, counts, num_genes)
    # Plot gene expression tragectory
    print('Generating Gene Expression Trajectory ...')
    plot_trajectory(out_dir, counts, num_genes)
    counts.to_csv('data/counts_clust.txt', index=False, sep='\t')

    print('Clustering Genes (this might take awhile) ...')
    subprocess.call(['mkdir', 'clust_out'])
    subprocess.call(['clust', 'data/counts_clust.txt', '-o', './clust_out'])


if __name__ == '__main__':
    main()
