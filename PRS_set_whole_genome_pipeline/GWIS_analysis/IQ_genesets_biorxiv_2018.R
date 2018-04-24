IQ_gene_sets <- read.csv("~/Documents/COGS_pathways_results/Gene_sets_to_be_used/IQ3_gene-sets.csv", stringsAsFactors = F)
IQ_gene_sets <- IQ_gene_sets[-c(1:3), c(1,2)]
colnames(IQ_gene_sets) <- IQ_gene_sets[1,]
IQ_gene_sets <- IQ_gene_sets[-1,]


write
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
library(data.table)
write.table(IQ_gene_sets[,2], file = "~/Documents/IQ_sets_to_be_converted.txt",quote = F,col.names = F,row.names = F)

entrez_gene_sets <- scan(file = "~/Documents/IQ_sets_converted.txt")
Entrez_gene_sets_gprof <- read.csv("~/Downloads/gprofiler_results_1029030649231.csv", stringsAsFactors = F)

ensembl=useDataset("hsapiens_gene_ensembl",mart = useMart("ensembl"))

data = IQ_gene_sets[,2]
ans <- unique(getBM(attributes = c("wikigene", "entrezgene"),    
                    filters = "wikigene",
                    values = data,
                    mart = ensembl) )