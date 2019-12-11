import argparse
import sys
import os
import pandas as pd
import pybedtools
import subprocess


def gene_fetcher(gene_list, promoter_output, bed_df,
                 upstream, downstream, fasta_file):
    # takes an input list and fetches the
    # genomic location, returning the promoter
    # of each gene

    # BED6 format
    gene_name_index = 3
    chr_index = 0
    start_index = 1
    stop_index = 2
    strand_index = 5

    genes = gene_list
    output = open(promoter_output, 'w')
    i = 0
    for line in genes:
        try:
            g = line.rstrip()
        except AttributeError:
            print(line)
        # First line only has one common entry called 'Genes'
        if i == 0:
            i += 1
            continue

        # Start from the second line idx=1
        try:
            gene_idx = bed_df.loc[bed_df['acc'] == g].index[0]
        except IndexError:
            print(g + ' not found')
            gene_idx = -1
            continue

        if gene_idx != -1:
            # if the gene matches, write the promoter to a file
            if bed_df.iloc[gene_idx, strand_index] == '-':
                prom_start = int(bed_df.iloc[gene_idx, stop_index]) + upstream
                prom_end = int(bed_df.iloc[gene_idx, stop_index]) - downstream
                if prom_end < 1:
                    prom_end = 1
                string = bed_df.iloc[gene_idx, chr_index] + '\t' +\
                    str(prom_end) + '\t' + str(prom_start) +\
                    '\t' + g + '\t' + '1' + '\t' + '-' + '\n'
                a = pybedtools.BedTool(string, from_string=True)
                a = a.sequence(fi=fasta_file)
                final = open(a.seqfn).read()
                output.write(final)
            if bed_df.iloc[gene_idx, strand_index] == '+':
                prom_start = int(bed_df.iloc[gene_idx, start_index]) - upstream
                if prom_start < 1:
                    prom_start = 1
                prom_end = int(bed_df.iloc[gene_idx, start_index]) + downstream
                string = bed_df.iloc[gene_idx, chr_index] + '\t' +\
                    str(prom_start) + '\t' + str(prom_end) +\
                    '\t' + g + '\t' + '1' + '\t' + '+' + '\n'
                a = pybedtools.BedTool(string, from_string=True)
                a = a.sequence(fi=fasta_file)
                final = open(a.seqfn).read()
                output.write(final)
        else:
            print('Gene: '+g+' not found.')
            continue
    output.close()
    # sys.exit(0)


def run_motif_enrichment(clust_out, ref_bed, ref_fa,
                         ref_motif, ame, up_dist=100, down_dist=100):
    '''
    Wraper function that parse out genes and
    extract promoter sequence from reference files
    '''
    df = pd.read_csv(clust_out, sep='\t')
    bed_df = pd.read_csv(ref_bed, sep='\t', header=None)
    bed_df = bed_df.iloc[:, 0:6]
    bed_df.columns = ['chr', 'start', 'end', 'acc', 'score', 'strand']
    bed_df['acc'] = bed_df['acc'].str.split('.', n=1, expand=True)[0]

    for i in df.columns:
        a = df[i].tolist()
        cleanedlist = [x for x in a if str(x) != 'NaN']
        cleanedlist = [x for x in cleanedlist if str(x) != 'nan']
        outname = str(i).split(' ')[0]
        gene_fetcher(cleanedlist, outname+'_promoters.bed',
                     bed_df, up_dist, down_dist, ref_fa)
        outfile = outname + '_promoters.bed'
        if ame:
            try:
                subprocess.check_call(['./ame_runner.sh', outfile, ref_motif])
            except ImportError:
                print('MEME Suite is not installed properly!')
        else:
            continue
    print('Analysis complete!')
