import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


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
    ax.scatter(finalDf.loc[:, 'PC1'],
               finalDf.loc[:, 'PC2'],
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
                            yticklabels=gene_names,
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
    '''
    Function to plot gene expression trajectory
    '''
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
