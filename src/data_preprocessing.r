if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install(version = "3.10")

# Load or install required packages
if (!requireNamespace("affy", quietly = TRUE)) {
    BiocManager::install("affy")
    require("affy")}
if (!requireNamespace("oligo", quietly = TRUE)) {
    BiocManager::install("oligo")
    require("oligo")}
if (!requireNamespace("limma", quietly = TRUE)) {
    BiocManager::install("limma")
    require("limma")}

if (!requireNamespace("mogene10sttranscriptcluster.db", quietly = TRUE)) {
    BiocManager::install("mogene10sttranscriptcluster.db")
    require("mogene10sttranscriptcluster.db")}

# function to extract mapped affyix matrix probe countent from db
extract <- function(x){
    as.list(x[mappedkeys(x)])}

setwd('./data')
celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles)

eset <- rma(affyRaw)

my_frame <- data.frame(exprs(eset))

Annot_accnum <- data.frame(sapply(extract(mogene10sttranscriptclusterACCNUM), paste, collapse=", "))
Annot_symbol <- data.frame(sapply(extract(mogene10sttranscriptclusterSYMBOL), paste, collapse=", "))

print('Annotating accession numbers ...')
all <- merge(Annot_accnum, my_frame, by.x="row.names", by.y="row.names", all.x=T)
print('Annotating gene symbols ...')
all <- merge(Annot_symbol, all, by.x="row.names", by.y=1, all.x=T)
colnames(all)[1:3] = c('AFFY_ID', 'SYMBOL','ACCNUM')

colnames(all)[4:(dim(all)[2])] = gsub(".*SP.","",colnames(all)[4:(dim(all)[2])])
colnames(all)[4:(dim(all)[2])] = gsub(".CEL","",colnames(all)[4:(dim(all)[2])])

print('Saving count table ...')
write.table(all,file="raw_counts.txt",sep="\t", row.names=F, col.names=T, quote=F)
