## calculating the similarity between the PGC data and CLOZUK data###
## Checking the removal of SNPs from the database

## prerequisites
library(data.table)

## setting working directory
wd<-getwd()

## Reading in PGC data
## Select for CHR 22
PGC_data <- fread(paste(wd,"/daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt",sep = ""))
test_PGC_data_chr_22_integer <- which(PGC_data[,CHR] == '22')
PGC_test_data_frame <- PGC_data[test_PGC_data_chr_22_integer,]


## Reading in example CLOZUK file
CLOZUK_data<- fread(paste(wd,"/CLOZUK_GWAS_BGE_chr22.bim",sep=""))
CLOZUK_test_data_frame<-CLOZUK_data[1:100,]

## checking overlap between both CLOZUK and schizophrenia data
PGC_integer_positions<- match(PGC_test_data_frame[,BP],CLOZUK_data[,V4])

test2<- which(test1 != 'NA')

CLOZUK_data[11,]
test3 <- which(PGC_test_data_frame[,BP] == '16051249')

