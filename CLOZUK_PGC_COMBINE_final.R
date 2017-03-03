#!/bin/env R
#################################################################
##################### PROPER SCRIPT START #######################
#################################################################

### Start Timer
ptm <- proc.time()


###############################
# Checking location for serial or batch analysis #
###############################

System_information <- Sys.info()
whereami <- System_information['user']

if (whereami == "johnhubert") {
  chromosome.number <- 22
  
} else if (whereami == 'JJ') {
  chromosome.number <- 22
  
} else if (whereami == "c1020109") {

# Preparing to run in Job array
AI <- Sys.getenv("PBS_ARRAY_INDEX")
chromosome.number <- as.numeric(AI)

} else {
  stop("current environment NOT at work/home or on servers, please add to script above to specify where you are and what type of analysis you want to do")
}


###############################
# ODDS RATIO TO BETA FUNCTION #
###############################
# add new environment#
e <- new.env()

# BETA and OR change function ###
log.to.odds <- function(imported.data.table) {
  imported.dt.col.names <- colnames(imported.data.table)
  if (any("OR" == imported.dt.col.names) == F) {
    cat("No Odds Ratio included in",deparse(substitute(imported.data.table)))
  }else{
    assign("PGC.BETA", log(imported.data.table$OR), envir = e)
  }
}

change.odds <- function (odds.ratios) {
  PGC.NEW.OR <- 1 / odds.ratios
  assign("PGC.NEW.OR", PGC.NEW.OR, envir = e)
}

##################################################################
#### CHECK TO SEE IF TWO COLUMNS HAVE ALLELES SWAPPED FUNCTION ###
##################################################################

Checking_allele_swapping <- function(alteredtable1,table2,which.is.combined = c("NONE","PGC","CLOZUK")){
  if (which.is.combined == "PGC") {
    Allele.x <- alteredtable1[,c("A1.x","A2.x"),with = F]
    Allele.y <- table2[,c("A1","A2"), with = F]
  }
  if(which.is.combined == "CLOZUK") {
    Allele.x <- alteredtable1[,c("A1.y","A2.y"),with = F]
    Allele.y <- table2[,c("A1","A2"), with = F]
  }
  if(which.is.combined == "NONE"){
    Allele.x <- alteredtable1[,c("A1","A2"),with = F]
    Allele.y <- table2[,c("A1","A2"), with = F]
  }
  
  Allele.x <- as.data.frame(Allele.x)
  Allele.y <- as.data.frame(Allele.y)
  iterations.to.remain <- NULL
  
  for (i in 1:length(Allele.x[,1])){
    if (Allele.x[i,1] != Allele.y[i,1]) {
      iterations.to.remain <- c(iterations.to.remain,i)
    }
  }
  argument.name <-deparse(substitute(alteredtable1))
  assign(paste0(argument.name,"flipped.alleles"), iterations.to.remain, envir = e)
}


#################################
# COMBINING CLOZUK AND PGC SNPS #
#################################

# load libraries
library(data.table)

# Set Working directory
setwd(".")
getwd()

### Adding in PGC data ###
PGC.data.frame <- fread(paste0("PGC_table",chromosome.number,".txt"))
PGC.data.frame.original <- copy(PGC.data.frame)

cat("Number of SNPS in PGC Chr:",chromosome.number, nrow(PGC.data.frame))

### Read in the CLOZUK data ###
CLOZUK.data <- fread(paste0("CLOZUK_GWAS_BGE_chr",chromosome.number,".bim"))
CLOZUK.original <- copy(CLOZUK.data)
# Number of SNPs in CLOZUK data
cat("Number of SNPs in CLOZUK Chr:",chromosome.number, nrow(CLOZUK.data))

### Replace the column names ###
setnames(CLOZUK.data, c("CHR","SNP","GENEDIST","BP","A1","A2"))

### Replacing PGC values###
PGC.integers.to.change <- PGC.data.frame[,.I[grep(paste0("^",chromosome.number,":\\d+:\\w+:\\w+"),SNP, perl = T,invert = T)]]
PGC.divisive.names <- PGC.data.frame[PGC.integers.to.change]
PGC.alternative <- PGC.data.frame[,SNP := paste0(CHR,":",BP)]
PGC.data.frame <- PGC.data.frame[,SNP := paste0(CHR,":",BP)]
PGC.alternative <- PGC.alternative[,c("CHR","SNP","BP","A1","A2","OR"),with = F]

### Adding BETA column to data
log.to.odds(PGC.alternative)
PGC.alternative[,BETA := e$PGC.BETA]
PGC.data.frame[,BETA := e$PGC.BETA]


### Changes the CLOZUK identifiers ###
CLOZUK.integers.to.change <- CLOZUK.data[,.I[grep(paste0("^",chromosome.number,":\\d+:\\w+:\\w+"),SNP, perl = T,invert = T)]]
CLOZUK.divisive.names <- CLOZUK.data[CLOZUK.integers.to.change]
CLOZUK.alternative <- CLOZUK.data[,SNP := paste0(CHR,":",BP)]
CLOZUK.alternative <- CLOZUK.alternative[,c("CHR","SNP","BP","A1","A2"),with = F]

# Simple tracking of each row for future checks
CLOZUK.alternative$place <- c(1:length(CLOZUK.alternative$CHR))
PGC.alternative$place <- c(1:length(PGC.alternative$CHR))

### Merging based on BP position ###
combined.CLOZUK.PGC <- merge(x = PGC.alternative,y = CLOZUK.alternative,by = c('BP','CHR'), all = F, sort = F)

# Change all alleles to uppercase: saving time method
# below it states:
## assign columns you want to change to object cols
## take all rows of data.table; apply the toupper function to Subset of Data columns named in cols object
cols <- c('A1.x','A2.x','A1.y','A2.y')
combined.CLOZUK.PGC[, (cols) := lapply(.SD, toupper), .SDcols = cols ]

cat("Number of SNPs BEFORE flipping between CLOZUK and PGC Chr:",chromosome.number ,nrow(combined.CLOZUK.PGC))

a <- which(combined.CLOZUK.PGC$A1.x == "C" & combined.CLOZUK.PGC$A2.x == "G"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.x == "G" & combined.CLOZUK.PGC$A2.x == "C"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.x == "A" & combined.CLOZUK.PGC$A2.x == "T"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.x == "T" & combined.CLOZUK.PGC$A2.x == "A"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}

cat("Chr:",chromosome.number ,nrow(combined.CLOZUK.PGC))

a <- which(combined.CLOZUK.PGC$A1.y=="C" & combined.CLOZUK.PGC$A2.y=="G"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.y=="G" & combined.CLOZUK.PGC$A2.y=="C"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.y=="A" & combined.CLOZUK.PGC$A2.y=="T"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.y=="T" & combined.CLOZUK.PGC$A2.y=="A"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}


cat("CLOZUK_PGC Chr:",chromosome.number, "remove A-T, C-G: N=", nrow(combined.CLOZUK.PGC))

## FLIPPING ##
# A1.x/A2.x are PGC alleles
# A1.y/A2.y are CLOZUK alleles

a <- which (combined.CLOZUK.PGC$A1.y == combined.CLOZUK.PGC$A2.x & combined.CLOZUK.PGC$A2.y == combined.CLOZUK.PGC$A1.x)
combined.CLOZUK.PGC$BETA[a] <- (-combined.CLOZUK.PGC$BETA[a])

# Change the values of the OR to suit Allele changes #
if (length(a) >= 1) {
  PGC.OR <- combined.CLOZUK.PGC$OR[a]
  change.odds(PGC.OR)
  combined.CLOZUK.PGC$OR[a]<- e$PGC.NEW.OR
  rm(PGC.NEW.OR, envir = e)
  rm(PGC.OR)
  combined.CLOZUK.PGC$A1.x[a] <- as.character(combined.CLOZUK.PGC$A1.y[a])
  combined.CLOZUK.PGC$A2.x[a] <- as.character(combined.CLOZUK.PGC$A2.y[a])
}

# Find alleles that are swapped in a more complicated manner (eg: with A -> G)
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

# Same OR change #
if (length(a) >= 1){
  combined.CLOZUK.PGC$BETA[a] <- (-combined.CLOZUK.PGC$BETA[a])
  PGC.OR <- combined.CLOZUK.PGC$OR[a]
  change.odds(PGC.OR)
  combined.CLOZUK.PGC$OR[a]<- e$PGC.NEW.OR
  rm(PGC.NEW.OR, envir = e)
  rm(PGC.OR)
  combined.CLOZUK.PGC$A1.x[a] <-as.character(combined.CLOZUK.PGC$A1.y[a])
  combined.CLOZUK.PGC$A2.x[a] <-as.character(combined.CLOZUK.PGC$A2.y[a])
}

a <- which(combined.CLOZUK.PGC$A1.y != combined.CLOZUK.PGC$A1.x | combined.CLOZUK.PGC$A2.y!=combined.CLOZUK.PGC$A2.x)
if (length(a)>0) combined.CLOZUK.PGC <- combined.CLOZUK.PGC[-a,]
cat("Chr:", chromosome.number,", remove A-T, C-G and other mismatches: N=", nrow(combined.CLOZUK.PGC))

### Find  out where in the original table the SNP identifiers are
PGC.final.index <- combined.CLOZUK.PGC$place.x
PGC.A1.allele.changes <- combined.CLOZUK.PGC$A1.x
PGC.A2.allele.changes <- combined.CLOZUK.PGC$A2.x
PGC.OR.changes <- combined.CLOZUK.PGC$OR
PGC.BETA.changes <- combined.CLOZUK.PGC$BETA
CLOZUK.final.index <- combined.CLOZUK.PGC$place.x


### Alter original files to have the right allele switching and the right BETA coefficient #
PGC.data.frame <- PGC.data.frame[PGC.final.index,A1:= PGC.A1.allele.changes]
PGC.data.frame <- PGC.data.frame[PGC.final.index,A2:= PGC.A2.allele.changes]
PGC.data.frame <- PGC.data.frame[PGC.final.index,OR:= PGC.OR.changes]
PGC.data.frame <- PGC.data.frame[PGC.final.index,BETA:= PGC.BETA.changes]


##### TODO #####
## Run in a job array using the job number as the chromosome number
## Change Valentina's script so that you are not changing the data within the data table...
## Find out what $B is and why it is set to negative...might be the negative strand for bp?
## Find Andrew's email about running on scratch and then run on scratch
## for the comments on $BP just change the effect size (stat in the PGC database) (changes the effect size based on the negative to positive)
## Change BETA, also change the table so that PGC are the alleles being changed rather than CLOZUK, because then you can change the effect sizes
## calculate BETA for PGC as well and incorporate that into this script
## Make the script more adept for other datasets as well

############################################
### CHECKING ABOVE ANALYSIS FOR ACCURACY ###
############################################

### checking everything is above board ###
Checking_allele_swapping(PGC.data.frame,PGC.data.frame.original,which.is.combined = "NONE")
Checking_allele_swapping(CLOZUK.alternative,CLOZUK.data,which.is.combined = "NONE")

if (length(e$CLOZUK.alternativeflipped.alleles) != 0) {
  warning("alleles have been flipped on CLOZUK")
}

checking.which.OR.do.not.equal.each.other <- which(PGC.data.frame$OR != PGC.data.frame.original$OR)
rm(PGC.data.frame.original)

problematic.SNPs <- match(e$PGC.data.frameflipped.alleles, checking.which.OR.do.not.equal.each.other)  
problematic.SNPs <- which(is.na(problematic.SNPs))
check1_index <- e$PGC.data.frameflipped.alleles[problematic.SNPs]
check1 <- PGC.data.frame[check1_index]
check2 <- PGC.data.frame[check1_index]

if (all(check1$OR == 1) & all(check2$OR == 1)){
  cat("You have",length(check1$CHR),"flipped SNPs with OR = 1")
} else {
  warning("There is an uneven amount of flipped SNPs")
}

### Checking conversions of SNPs ###
checking.duplications.PGC <- which (duplicated(combined.CLOZUK.PGC$SNP.x))
checking.duplications.CZK <- which (duplicated(combined.CLOZUK.PGC$SNP.y))

if (length(checking.duplications.CZK) != 0 & length(checking.duplications.PGC) != 0) {
  warning("There are duplicated SNPs common between CLOZUK and PGC")
}

## Altering data tables to fit previous input files ##
<<<<<<< HEAD
PGC.data.frame <- PGC.data.frame[, OR := NULL]
=======
PGC.data.frame <- PGC.data.frame[,-"OR",with = F]
>>>>>>> f1453e0ee59bc2ed7d1165fcb816ea27845efcb2
namesPGC <- names(PGC.data.frame)
namesPGC <- namesPGC[c(1:8,17,9:16)]
PGC.data.frame <- PGC.data.frame[, namesPGC, with = F]

## checking for combined SNPs ##
SNPs <- match(combined.CLOZUK.PGC$SNP.x,combined.CLOZUK.PGC$SNP.y)
if (length(which(is.na(SNPs))) > 0 ) {
  warning("there is an uneven number of SNPs common between both data.tables")
}

# Write out update file for use in plink
CLOZUK_together <- CLOZUK.original[,.(V2)][,SNP := CLOZUK.alternative$SNP]
rm(CLOZUK.original)

# Write according to destination
if (whereami == 'johnhubert' | whereami == 'JJ'){
  filename.CLOZUK.together <- "./output/CLOZUK_chr22_chr.pos.txt"
  new.PGC.table <- "./output/PGC_table22_new.txt"
  filename.common.snps <- "./output/chr22PGC_CLOZUK_common_SNPs.txt"
} else {  
  filename.CLOZUK.together <- paste0("./output/CLOZUK_chr",chromosome.number,"_chr.pos.txt")
  new.PGC.table <- paste0("./output/PGC_table",chromosome.number,"_new.txt")
  filename.common.snps <- paste0("./output/chr",chromosome.number,"PGC_CLOZUK_common_SNPs.txt")
}

# Write update file for plink
write.table(CLOZUK_together, file = filename.CLOZUK.together, quote = F, row.names = F, col.names = F)

# Write out new PGC table
write.table(PGC.data.frame, file = new.PGC.table, row.names = F, quote = F)

## Write out the common SNPs between CLOZUK and PGC
write(combined.CLOZUK.PGC$SNP.x,file = filename.common.snps)

### End timer
proc.time() - ptm


#################################################################
##################### PROPER SCRIPT END #########################
#################################################################

### TODO ###

## Incorporate the new ID's into CLOZUK data including the genotype information etc... tar all together and keep in bim format#
## correspond SNP's to genomic regions by reading in the NCBI data source

#---------------------------
## clumping per gene ##
## plink --bfile CLOZUK_GWAS_BGE_chr22 --clump CLOZUK2_COGS_GWAS_noPGC2.assoc.dosage --clump-r2 0.1 --clump-kb 3000 --out test3 --clump-p1 0.0001 --clump-verbose --clump-range input MAGMA limits
## rewrite all of this, use the altered.names files to get the right list of SNPs for the CLOZUK dataset, PGC is fine, you cacn move it around alright


