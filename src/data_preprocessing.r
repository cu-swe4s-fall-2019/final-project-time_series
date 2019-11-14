source("http://bioconductor.org/biocLite.R")


# Load or install required packages
if (require("affy") == FALSE) {
    BiocInstaller::biocLite("affy")
    require("affy")}
if (require("oligo") == FALSE) {
    BiocInstaller::biocLite("oligo")
    require("oligo")}
if (require("limma") == FALSE) {
    BiocInstaller::biocLite("limma")
    require("limma")}

if (require("mogene10sttranscriptcluster.db") == FALSE) {
    BiocInstaller::biocLite("mogene10sttranscriptcluster.db")
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
