#!/bin/bash
mkdir ./data
echo Downloading raw Affymetrix data file ...
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75195/suppl/GSE75195_RAW.tar
tar -xvf GSE75195_RAW.tar -C ./data; rm GSE75195_RAW.tar
gunzip ./data/*.gz

# one of the files is missing a period separation
mv ./data/GSM1944906_EA07068_332444_MOGENE-1_0-ST-V1_12.BB.SP3H.IFN_1.CEL \
./data/GSM1944906_EA07068_332444_MOGENE.1_0.ST.V1_12.BB.SP.3H.IFN_1.CEL

echo Converting and annotating the count table ...
rscript ./src/data_preprocessing.r
