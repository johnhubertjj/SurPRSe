#################################################################
##################### PROPER SCRIPT START #######################
#################################################################


#################################
# COMBINING CLOZUK AND PGC SNPS #
#################################

# Preparing to run in Job array
AI <- Sys.getenv("PBS_ARRAY_INDEX")
chromosome.number <- as.numeric(AI)

# load libraries
library(data.table)

# Set Working directory
setwd(".")

### Adding in PGC data ###
PGC.test.data.frame <- fread(paste0("gzip -dc PGC_table",chromosome.number,".txt.gz"))

### Read in the CLOZUK data ###
CLOZUK.data <- fread(paste0("CLOZUK_GWAS_BGE_chr",chromosome.number,".bim"))

### Replace the column names ###
setnames(CLOZUK.data, c("CHR","SNP","Something","BP","A1","A2"))

###Replacing PGC values###
PGC.test.data.frame2 <- PGC.test.data.frame
PGC.integers.to.change <- PGC.test.data.frame2[,.I[grep(paste0(chromosome.number,":\\d+:\\w+:\\w+"),SNP, perl = T,invert = T)]]
PGC.divisive.names <- PGC.test.data.frame2[PGC.integers.to.change]
PGC.alternative <- PGC.test.data.frame2[PGC.integers.to.change, SNP := paste0(CHR,":",BP,":",A2,":",A1)]
PGC.alternative <- PGC.alternative[,c("CHR","SNP","BP","A1","A2"),with = F]
rm(PGC.test.data.frame2)


### Changes the CLOZUK identifiers in 3 lines of a Data.table in no time whatsoever ###
CLOZUK.data2 <- CLOZUK.data
CLOZUK.integers.to.change <- CLOZUK.data2[,.I[grep(paste0("^",chromosome.number,":\\d+:\\w+:\\w+"),SNP, perl = T,invert = T)]]
CLOZUK.divisive.names <- CLOZUK.data2[CLOZUK.integers.to.change]
CLOZUK.alternative <- CLOZUK.data2[CLOZUK.integers.to.change, SNP := paste0(CHR,":",BP,":",A2,":",A1)]
CLOZUK.alternative <- CLOZUK.alternative[,c("CHR","SNP","BP","A1","A2"),with = F]
rm(CLOZUK.data2)

### Merging based on BP position ###
combined.CLOZUK.PGC <- merge(CLOZUK.alternative,PGC.alternative,by = c('BP','CHR'), all = F, sort = F)

combined.CLOZUK.PGC$A1.x <- toupper(combined.CLOZUK.PGC$A1.x)
combined.CLOZUK.PGC$A2.x <- toupper(combined.CLOZUK.PGC$A2.x)
combined.CLOZUK.PGC$A1.y <- toupper(combined.CLOZUK.PGC$A1.y)
combined.CLOZUK.PGC$A2.y <- toupper(combined.CLOZUK.PGC$A2.y)

cat("Chr:",chromosome.number ,nrow(combined.CLOZUK.PGC))

a <- which(combined.CLOZUK.PGC$A1.x == "C" & combined.CLOZUK.PGC$A2.x == "G"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.x == "G" & combined.CLOZUK.PGC$A2.x == "C"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.x == "A" & combined.CLOZUK.PGC$A2.x == "T"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.x == "T" & combined.CLOZUK.PGC$A2.x == "A"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}

cat("Chr:",chromosome.number ,nrow(combined.CLOZUK.PGC))

a <- which(combined.CLOZUK.PGC$A1.y=="C" & combined.CLOZUK.PGC$A2.y=="G"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.y=="G" & combined.CLOZUK.PGC$A2.y=="C"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.y=="A" & combined.CLOZUK.PGC$A2.y=="T"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.y=="T" & combined.CLOZUK.PGC$A2.y=="A"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}


cat("Chr:",chromosome.number, "remove A-T, C-G: N=", nrow(combined.CLOZUK.PGC))

## FLIPPING ##
# A1.x/A2.x are CLOZUK alleles
# A1.y/A2.y are PGC alleles

a <- which (combined.CLOZUK.PGC$A1.y == combined.CLOZUK.PGC$A2.x & combined.CLOZUK.PGC$A2.y == combined.CLOZUK.PGC$A1.x)
# combined.CLOZUK.PGC$BP[a] <- (-combined.CLOZUK.PGC$BP[a])
combined.CLOZUK.PGC$A1.x[a] <- as.character(combined.CLOZUK.PGC$A1.y[a])
combined.CLOZUK.PGC$A2.x[a] <- as.character(combined.CLOZUK.PGC$A2.y[a])

a <- which(combined.CLOZUK.PGC$A1.y == combined.CLOZUK.PGC$A1.x & combined.CLOZUK.PGC$A2.y == combined.CLOZUK.PGC$A2.x)
d <- seq(1:nrow(combined.CLOZUK.PGC));d <- d[-a]

b1 <- which(combined.CLOZUK.PGC$A1.x[d] == "A")
b2 <- which(combined.CLOZUK.PGC$A1.x[d] == "C")
b3 <- which(combined.CLOZUK.PGC$A1.x[d] == "G")
b4 <- which(combined.CLOZUK.PGC$A1.x[d] == "T")

combined.CLOZUK.PGC$A1.x[d[b1]] <- "T"
combined.CLOZUK.PGC$A1.x[d[b2]] <- "G"
combined.CLOZUK.PGC$A1.x[d[b3]] <- "C"
combined.CLOZUK.PGC$A1.x[d[b4]] <- "A"

b1 <- which(combined.CLOZUK.PGC$A2.x[d] == "A")
b2 <- which(combined.CLOZUK.PGC$A2.x[d] == "C")
b3 <- which(combined.CLOZUK.PGC$A2.x[d] == "G")
b4 <- which(combined.CLOZUK.PGC$A2.x[d] == "T")

combined.CLOZUK.PGC$A2.x[d[b1]] <-"T"
combined.CLOZUK.PGC$A2.x[d[b2]] <-"G"
combined.CLOZUK.PGC$A2.x[d[b3]] <-"C"
combined.CLOZUK.PGC$A2.x[d[b4]] <-"A"

a <- which(combined.CLOZUK.PGC$A1.y==combined.CLOZUK.PGC$A2.x & combined.CLOZUK.PGC$A2.y==combined.CLOZUK.PGC$A1.x)
# combined.CLOZUK.PGC$BP[a] <- (-combined.CLOZUK.PGC$BP[a])
combined.CLOZUK.PGC$A1.x[a] <-as.character(combined.CLOZUK.PGC$A1.y[a])
combined.CLOZUK.PGC$A2.x[a] <-as.character(combined.CLOZUK.PGC$A2.y[a])

a <- which(combined.CLOZUK.PGC$A1.y != combined.CLOZUK.PGC$A1.x | combined.CLOZUK.PGC$A2.y!=combined.CLOZUK.PGC$A2.x)
if (length(a)>0) combined.CLOZUK.PGC <- combined.CLOZUK.PGC[-a,]

cat("Chr:", chromosome.number,", remove A-T, C-G and other mismatches: N=", nrow(combined.CLOZUK.PGC))

z <- gzfile(paste0("PGC_CLOZUK_SNP_table",chromosome.number,".txt.gz"))
write.table(combined.CLOZUK.PGC,file = z, row.names = F, quote = F)


##### TODO #####

# Run in a job array using the job number as the chromosome number
# Change Valentina's script so that you are not changing the data within the data table
# Write shell script that can extract the bim files from the CLOZUK BGE dataset
# Find out what $B is and why it is set to negative...might be the negative strand for bp?
# 

