IQ_gene_sets <- read.csv("~/Documents/COGS_pathways_results/Gene_sets_to_be_used/IQ3_gene-sets.csv", stringsAsFactors = F)
IQ_gene_sets <- IQ_gene_sets[-c(1:3), c(1,2)]
colnames(IQ_gene_sets) <- IQ_gene_sets[1,]
IQ_gene_sets <- IQ_gene_sets[-1,]

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
library(data.table)
library(stringr)
write.table(IQ_gene_sets[,2], file = "~/Documents/IQ_sets_to_be_converted.txt",quote = F,col.names = F,row.names = F)

entrez_gene_sets <- scan(file = "~/Documents/IQ_sets_converted.txt")
Entrez_gene_sets_gprof <- read.csv("~/Downloads/gprofiler_results_1029030649231.csv", stringsAsFactors = F, header = F)

ensembl=useDataset("hsapiens_gene_ensembl",mart = useMart("ensembl"))

data = IQ_gene_sets[,2]
ans <- unique(getBM(attributes = c("wikigene", "entrezgene"),    
                    filters = "wikigene",
                    values = data,
                    mart = ensembl) )

Entrez_gene_sets_gprof <- read.csv("~/Downloads/gprofiler_results_1029030649231.csv", stringsAsFactors = F, header = F)
str(Entrez_gene_sets_gprof)
Entrez_gene_sets_gprof <- as.data.table(Entrez_gene_sets_gprof)
Entrez_gene_sets_gprof <- Entrez_gene_sets_gprof[,V4 := str_replace(Entrez_gene_sets_gprof$V4, "ENTREZGENE_ACC:","")]
Entrez_gene_sets_gprof <- Entrez_gene_sets_gprof[,V4 := str_replace(Entrez_gene_sets_gprof$V4, "\\bN/A\\b", "NA")]

setnames(Entrez_gene_sets_gprof, old = "V5", new = "Gene")
length(which(Entrez_gene_sets_gprof$V4 == "NA"))

setnames(Entrez_gene_sets_gprof, old = "V4", new = "entrez_gene_sets")
names(entrez_gene_sets) <- "entrez_gene_sets"
entrez_gene_sets <- as.data.frame(entrez_gene_sets)
Entrez_gene_sets_gprof$entrez_gene_sets <- as.numeric(Entrez_gene_sets_gprof$entrez_gene_sets)

testing_ID_s <- merge(Entrez_gene_sets_gprof, entrez_gene_sets, by = "entrez_gene_sets")

all_gene_sets <- unique(testing_ID_s$entrez_gene_sets)
testing_ID_s <- testing_ID_s[unique(entrez_gene_sets)]

names(all_gene_sets) <- "entrez_gene_sets"
all_gene_sets <- as.data.frame(all_gene_sets)

setnames(testing_ID_s, old = "V2", new = "Gene")

Final_gene_sets <- merge(IQ_gene_sets, testing_ID_s, by = "Gene")


# Enter the GMT style scripts

GENES_to_snps <- scan("~/Documents/COGS_cognition_paper/Data/Raw_data/c5.all.v6.1.entrez.gmt", what = "", sep = "\n")

## Convert from MAGMA structure to an object easily analysed in R ##
y <- strsplit(GENES_to_snps, "[[:space:]]+")
names(y) <- sapply(y, '[[', 1)
y <- lapply(y, '[', -1)
y <- lapply(y,'[', -1)

Melted_MAGMA_list <- melt(y)
names(Melted_MAGMA_list) <- c("Entrez_gene", "Gene-set")

Melted_MAGMA_list <- as.data.table(Melted_MAGMA_list)
Names_of_gene_sets <- unique(IQ_gene_sets$`Gene-set`)

setkey(Melted_MAGMA_list,"Gene-set")
Final_entrez_gene_set <- Melted_MAGMA_list[Names_of_gene_sets]
Final_entrez_gene_set <- merge(IQ_gene_sets,Melted_MAGMA_list,by ="Gene-set", all = F) 

setcolorder(Final_entrez_gene_set,neworder = c("Gene-set","Entrez_gene"))

write.table(Final_entrez_gene_set,"~/Documents/COGS_cognition_paper/Data/Raw_data/IQ_2018_biorxiv_genesets.txt",quote = F,col.names = F,row.names = F)

