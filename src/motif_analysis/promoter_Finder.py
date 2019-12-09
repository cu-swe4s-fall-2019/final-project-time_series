import argparse
import sys
import os
import pandas as pd
import pybedtools
import subprocess

def gene_fetcher(gene_list, promoter_output, bedfile, upstream, downstream, fasta_file):
    # takes an input list and fetches the genomic location, returning the promoter
    # of each gene

    genes = gene_list
    output = open(promoter_output, 'w')
    i = 0
    for line in genes:
        g = line.rstrip()
        if i == 0:
            i += 1
            continue
        f = open(bedfile)
        for l in f:
            # Standard BED format for bed6
            sample_line = l.rstrip().split('\t')
            gene_name_index = 3
            chr_index = 0
            start_index = 1
            stop_index = 2
            strand_index = 5
            gene = sample_line[gene_name_index]
            if gene == g:
                # if the gene matches, write the promoter to a file
                if sample_line[strand_index] == '-':
                    prom_start = int(sample_line[stop_index]) + upstream
                    prom_end = int(sample_line[stop_index]) - downstream
                    if prom_end < 1:
                        prom_end = 1
                    string =  sample_line[chr_index] + '\t' + str(prom_end) +\
                              '\t' + str(prom_start) + '\t' + gene + '\t' + '1' + '\t' + '-' + '\n'                 
                    a = pybedtools.BedTool(string, from_string=True)
                    a = a.sequence(fi = fasta_file)
                    final = open(a.seqfn).read()
                    output.write(final)
                if sample_line[strand_index] == '+':
                    prom_start = int(sample_line[start_index]) - upstream
                    if prom_start < 1:
                        prom_start = 1
                    prom_end = int(sample_line[start_index]) + downstream
                    string = sample_line[chr_index] +'\t' + str(prom_start) + '\t' + str(prom_end) +\
                             '\t' + gene + '\t' + '1' + '\t' + '+' + '\n'
                    a = pybedtools.BedTool(string, from_string=True)
                    a = a.sequence(fi = fasta_file)
                    final = open(a.seqfn).read()
                    output.write(final)
    output.close()
    # sys.exit(0)

def main():
    parser = argparse.ArgumentParser(
            description='fetch promoter information for input genes')

    parser.add_argument('--bedfile',
                        type=str,
                        help='path to BED file containing gene coordinates',
                        required=True)

    parser.add_argument('--upstream_distance',
                        type=int,
                        help='number of bases upstream to annotated start site\
                             for promoter determination',
                        default=100,
                        required=False)

    parser.add_argument('--downstream_distance',
                        type=int,
                        help='number of bases downstream to annotated start site\
                             for promoter determination',
                        default=100,
                        required=False)
    
    parser.add_argument('--input_genes',
                        type=str,
                        help='list of input genes',
                        required=True)

    parser.add_argument('--fasta',
                        type=str,
                        help='genome in fasta format',
                        required=True)
    parser.add_argument('--motif_file',
                        type=str,
                        help='motifs in MEME format',
                        required=True)

    args = parser.parse_args()
    
    
    df = pd.read_csv(args.input_genes, sep='\t')
    for i in df.columns:
        a = df[i].tolist()
        cleanedlist = [x for x in a if str(x) != 'nan']
        outname = str(i).replace(" ", "_")
        gene_fetcher(cleanedlist, outname + '_promoters.bed', args.bedfile, args.upstream_distance, args.downstream_distance, args.fasta)
        outfile = outname + '_promoters.bed'
        subprocess.check_call(['./ame_runner.sh', outfile, args.motif_file])


if __name__ == '__main__':
    main()
