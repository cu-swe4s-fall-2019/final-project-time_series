# Time Series Expression Analysis
This is a implementation of analysis that is geared towards analyzing time series expression data. The program will take raw read counts from time series expression data, normalize the data, cluster the genes based on expression patter, extract promoter sequences of each gene cluster and conduct motif enrichment analysis to figure out what underlying transcription factor is potentially regulating such gene expression.

## Usage
### Input count data
The data is downloaded from NCBI GEO database (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75306). The data was originally generated for the publication by [Mostafavi, S., et al.](https://doi.org/10.1016/j.cell.2015.12.032). To simplify the process, the raw microarray data is processed and normalized for you. The unprocessed microarray files along with the final raw counts output can be found in `./data/`.

### Input reference data
In order to extract necessary promoter sequences from a list of clustered genes, a series of reference files are needed. This includes:
- mm10.refseq.bed
- mm10.fasta

While the reference bed file is included within `./ref/`, the reference genome file is not due to its large size. So please download the reference genome for mouse:
```
wget -O ./ref/mm10.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.p6.genome.fa.gz
gunzip ./ref/mm10.fa.gz
```

### Data Processing and clustering
To process the raw read count data, a python script is included with default parameter. Once the input data is downloaded. Simply run the python script `main.py`:
```
python main.py
```
Additionally, this program is written in a way to handle any time-series expression data with different time points within each column and genes on each row. To customize specific parameters for your need, you can adjust the following parameters:

Arguments |  Description
--- | ---
-h  | Display help message
-o  | Specify the out directory to store plots <Default: ./out>
-c  | Input read counts file name <Default: ./data/raw_counts.txt>
-k  | Manually conduct K means clustering with specified K. This disables optimization. <Default: None>
-n  | Number of genes to include for generating heatmap and trajectory plots <Default: 50 >
-a  | Enable Analysis of motif enrichment (AME) <Default: False>

## Data
Mouse B cells were treated with Interferon alpha over a time course, the gene expression over time is normalized and visualized in the plot below:
![](./out/gene_heatmap.png)
![](./out/gene_trajectory.png)

## Clustering
K-means clustering optimization
![](./out/SSvK.png)
Clustering results
![](./out/cluster.png)

## Motif Enrichment
Cluster 1 gene promoter motif enrichment
![](./out/cluster1_enrichment.png)

Cluster 4 gene promoter motif enrichment
![](./out/cluster4_enrichment.png)
