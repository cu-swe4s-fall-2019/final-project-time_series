# Time Series Expression Analysis

## Data Downloading and Preprocessing
To download and process data, execute the shell script:
```
bash run.sh
```
This will download the raw data from NCBI GEO database, extract raw data from CEL files then summarize and annotate the files. The final raw counts output can be found in `./data/raw_counts.txt`

Download the reference genome for mouse:
```
wget -O ./ref/mm10.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.p6.genome.fa.gz
gunzip ./ref/mm10.fa.gz
```


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
