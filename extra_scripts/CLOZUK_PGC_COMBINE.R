
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
untar(paste0("CLOZUK_GWAS_BGE_chr",chromosome.number,".tar.gz"),files = paste0("CLOZUK_GWAS_BGE_chr",chromosome.number,".bim"))
CLOZUK.data <- fread(paste0("CLOZUK_GWAS_BGE_chr",chromosome.number,".bim"))

### Replace the column names ###
setnames(CLOZUK.data, c("CHR","SNP","STAT","BP","A1","A2"))

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

CLOZUK.alternative$place <- c(1:length(CLOZUK.alternative$CHR))
PGC.alternative$place <- c(1:length(PGC.alternative$CHR))
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
combined.CLOZUK.PGC$Effect[a] <- (-combined.CLOZUK.PGC$Effect[a])
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
combined.CLOZUK.PGC$Effect[a] <- (-combined.CLOZUK.PGC$Effect[a])
combined.CLOZUK.PGC$A1.x[a] <-as.character(combined.CLOZUK.PGC$A1.y[a])
combined.CLOZUK.PGC$A2.x[a] <-as.character(combined.CLOZUK.PGC$A2.y[a])
a <- which(combined.CLOZUK.PGC$A1.y != combined.CLOZUK.PGC$A1.x | combined.CLOZUK.PGC$A2.y!=combined.CLOZUK.PGC$A2.x)
if (length(a)>0) combined.CLOZUK.PGC <- combined.CLOZUK.PGC[-a,]

cat("Chr:", chromosome.number,", remove A-T, C-G and other mismatches: N=", nrow(combined.CLOZUK.PGC))

# Find  out where in the original table the SNP identifiers are
CLOZUK.final.index <- combined.CLOZUK.PGC$place.x
PGC.final.index <- combined.CLOZUK.PGC$place.y

# write the combined tables and the tables with altered SNP identifiers to the wd
za <- gzfile(paste0("PGC_CLOZUK_SNP_table",chromosome.number,".txt.gz"))
zb <- gzfile(paste0("CLOZUK_altered_names_chr",chromosome.number,".txt.gz"))
zc <- gzfile(paste0("PGC_altered_names_chr",chromosome.number,".txt.gz"))

write.table(combined.CLOZUK.PGC,file = za, row.names = F, quote = F)
write.table(CLOZUK.alternative,file = zb, row.names = F, quote = F)
write.table(PGC.alternative,file = zc,row.names = F, quote = F)

# write the index as well just in case for clumping
write(x = CLOZUK.final.index,file = paste0("CLOZUK.merged.index.chr",chromosome.number,".txt"))
write(x = PGC.final.index,file = paste0("PGC.merged.index.chr",chromosome.number,".txt"))


##### TODO #####
## Run in a job array using the job number as the chromosome number
## Change Valentina's script so that you are not changing the data within the data table...
## Find out what $B is and why it is set to negative...might be the negative strand for bp?
## Find Andrew's email about running on scratch and then run on scratch
## for the comments on $BP just change the effect size (stat in the PGC database) (changes the effect size based on the negative to positive)

