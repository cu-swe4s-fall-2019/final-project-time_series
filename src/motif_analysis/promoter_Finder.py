import argparse
import sys
import os


def gene_fetcher(gene_list, bedfile, upstream, downstream):
    # takes an input list and fetches the genomic location, returning the promoter
    # of each gene

    genes = open(gene_list)
    output = open("promoter_file.bed", 'w')
    for line in genes:
        g = line.rstrip()
        f = open(bedfile)
        for l in f:
            sample_line = l.rstrip().split('\t')
            gene_name_index = 0
            chr_index = 1
            start_index = 2
            stop_index = 3
            strand_index = 4
            gene = sample_line[gene_name_index]
            if gene == g:
                # if the gene matches, write the promoter to a file
                if sample_line[strand_index] == '-':
                    prom_start = int(sample_line[stop_index]) + upstream
                    prom_end = int(sample_line[stop_index]) - downstream
                    if prom_end < 1:
                        prom_end = 1
                    string = gene + '\t' + str(prom_end) + '\t' + str(prom_start) + '\t' + '-' + '\n'
                    output.write(string)
                if sample_line[strand_index] == '+':
                    prom_start = int(sample_line[start_index]) - upstream
                    if prom_start < 1:
                        prom_start = 1
                    prom_end = int(sample_line[start_index]) + downstream
                    string = gene + '\t' + str(prom_start) + '\t' + str(prom_end) + '\t' + '+' + '\n'
                    output.write(string)
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

    args = parser.parse_args()
    
    

    gene_fetcher(args.input_genes, args.bedfile, args.upstream_distance, args.downstream_distance)

if __name__ == '__main__':
    main()
