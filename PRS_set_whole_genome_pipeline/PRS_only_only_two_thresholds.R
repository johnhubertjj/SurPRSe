library(data.table)

##Read in PGC.data

PGC <- fread("./combined_BIPvsSCZ_with_CHR.POS_identifiers.txt")
test_data_frame <- fread("./CLOZUK_BIPvsSCZ_FULL_GENOME_CLUMPED.bim")
setnames(test_data_frame, c("CHR", "SNP", "BP", "Gene_ID", "BP_1", "BP_2"))
a <- merge(PGC,test_data_frame,by="SNP", all=F)
## Excess info if required ##
#############################
## Plink --gene-report tmp1.txt Chromosom_MAGMA_GENE_REGIONS_NCBI37.3.txt
##       --out CLOZUK_BGE_GENE_REGIONS

# Plink --bfile CLOZUK_GWAS_BGE_chr22
#       --make-set Chromosome_22_MAGMA_GENE_REGIONS_NCBI37.3.txt
#       --out testlol
#       --write-set
#############################

# Ideas for getting scores and keeping them within the R environment

# 1 use MAGMA and create the gene.annot file: read as a list into R 
# then create a Data.table with repeated SNP's for each gene in the format SNP, A1, BETA, GENE
# Then parse through each row using the data table notation to create a score and profile wihtin an index
# Problem here will be the output to plink
# Can maybe put into a list of data tables for the scoring? Find out if it is possible...if not then just input, print and delete
# 

# Scoring per gene: rough and undescriptive at the moment

# SCORING 
a1 <- a[,.(SNP,A1,BETA)]

a0.05 <-a [P <= 0.05, .(SNP,A1,BETA)]

filenamea1 <- "CLOZUK_BIPvsSCZ_FULL_GENOME_CLUMPED_1_threshold.score"
filenamea0.05 <- "CLOZUK_BIPvsSCZ_FULL_GENOME_CLUMPED_0.05_threshold.score"
write.table(file = filenamea1, a1, row.names = F, col.names = F, quote = F, sep="\t")
write.table(file = filenamea0.05, a0.05, row.names = F, col.names = F, quote = F, sep="\t")
