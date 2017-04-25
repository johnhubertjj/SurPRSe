#!/bin/env R
#################################################################
##################### PROPER SCRIPT START #######################
#################################################################

### Start Timer
ptm <- proc.time()

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)

##################################################
# Checking location for serial or batch analysis #
##################################################

System_information <- Sys.info()
whereami <- System_information['user']

if (whereami == "johnhubert") {
  chromosome.number <- args[7]
  
} else if (whereami == 'JJ') {
  chromosome.number <- args[7]
  
} else if (whereami == "c1020109") {

# Preparing to run in Job array
AI <- Sys.getenv("PBS_ARRAY_INDEX")
chromosome.number <- as.numeric(AI)

} else {
  stop("current environment NOT at work/home or on servers, please add to script above to specify where you are and what type of analysis you want to do")
}



# specify the different input tables #
Training_datatable <- paste0("./output/",args[3],"_new.txt")
Validation_datatable_bim_file <- paste0(args[4],".bim")
Training_name <- args[5]
Validation_name <- args[6]
chromosome.number <- args[7]
CHR <- args[8]
SNP <- args[9]
BP <- args[10]
A1 <- args[11]
A2 <- args[12]
OR <- args[13]
BETA <- args[14]
number_of_frequency_columns <- args[15]

Column_headers_to_keep <- c(CHR, SNP, BP, A1, A2, OR, BETA)
Column_headers <- c("CHR", "SNP", "BP", "A1", "A2", "OR", "BETA")
Column_headers_to_keep_elements <- which(Column_headers_to_keep == TRUE)
Column_headers <- Column_headers[Column_headers_to_keep_elements]

 #setwd("~/Documents/testing_cross_disorder") 
 #Validation_datatable_bim_file <- "CLOZUK_GWAS_BGE_chr14.bim"
# Validation_datatable_bim_file <- "/Volumes/HD-PCU2/g1000_European_cohort_only_NO_bad_LD_Extracted.bim"
# Training_datatable <- "/Volumes/HD-PCU2/ADD/daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta"
 #Training_datatable <- "/Volumes/HD-PCU2/PGC_noCLOZUK_data/PGC_table14.txt"
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

Checking_length_of_alleles <- function(input_table1){
  # Removes all multiple allele counts (aka indels and deletions) to clean up the data.
  input_table1 <- input_table1[!(nchar(A1.x) > 1 | nchar(A2.x) > 1)]
  input_table1 <- input_table1[!(nchar(A1.y) > 1 | nchar(A2.y) > 1)]
  
  # Checks for any Alleles which are not standard 1:1
  integers_to_check <- which(input_table1$A1.x != "A" & input_table1$A1.x != "C" & input_table1$A1.x !=  "T" & input_table1$A1.x !=  "G")
  integers_to_check2 <- which(input_table1$A1.y != "A" & input_table1$A1.y != "C" & input_table1$A1.y !=  "T" & input_table1$A1.y !=  "G")
  integers_to_check3 <- which(input_table1$A2.x != "A" & input_table1$A2.x != "C" & input_table1$A2.x !=  "T" & input_table1$A2.x !=  "G")
  integers_to_check4 <- which(input_table1$A2.y != "A" & input_table1$A2.y != "C" & input_table1$A2.y !=  "T" & input_table1$A2.y !=  "G")
   
  # Would probably work better in a loop but I got too lazy
  if(length(integers_to_check) != 0){
    input_table1 <- input_table1[!integers_to_check]
  }
  
  if(length(integers_to_check2) != 0){
    input_table1 <- input_table1[!integers_to_check2]
  }
  
  if(length(integers_to_check3) != 0){
    input_table1 <- input_table1[!integers_to_check3]
  }
  
  if(length(integers_to_check4) != 0){
    input_table1 <- input_table1[!integers_to_check4]
  }
  assign("combined.CLOZUK.PGC", input_table1, envir = .GlobalEnv)
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
PGC.data.frame <- fread(Training_datatable)
PGC.data.frame.original <- copy(PGC.data.frame)
original_column_headers_PGC <- colnames(PGC.data.frame.original)

cat("Number of SNPS in ",Training_name," Chr:",chromosome.number, "N=" ,nrow(PGC.data.frame))

### Remove Duplicated SNPs in Training here ###
PGC.duplicated.removed <- which(duplicated(PGC.data.frame$SNP))
PGC.duplicated.removed.rev <- which(duplicated(PGC.data.frame$SNP, fromLast = T))

### One-line duplication check - Training ###
length(PGC.data.frame$SNP)
if (length(PGC.duplicated.removed) >= 1){
  PGC_duplicate_SNPS <- PGC.data.frame$SNP[(duplicated(PGC.data.frame$SNP) | duplicated(PGC.data.frame$SNP, fromLast = TRUE)) ]
  PGC.data.frame <- PGC.data.frame[!(duplicated(PGC.data.frame$SNP) | duplicated(PGC.data.frame$SNP, fromLast = TRUE)) ]
}
length(PGC.data.frame$SNP)

cat("Number of SNPS in ",Training_name," after duplication check Chr:",chromosome.number, "N=" ,nrow(PGC.data.frame))

### Read in the CLOZUK data ###
CLOZUK.data <- fread(Validation_datatable_bim_file)
CLOZUK.original <- copy(CLOZUK.data)
# Number of SNPs in CLOZUK data
cat("Number of SNPs in", Validation_name, "Chr:",chromosome.number, "N=" ,nrow(CLOZUK.data))

### Replace the column names ###
setnames(CLOZUK.data, c("CHR","SNP","GENEDIST","BP","A1","A2"))

### Remove Duplicated SNPs in Validation here ###
CLOZUK.duplicated.removed <- which(duplicated(CLOZUK.data$SNP))
CLOZUK.duplicated.removed.rev <- which(duplicated(CLOZUK.data$SNP, fromLast = T))

### One-line duplication check - Training ###
length(CLOZUK.data$SNP)
if (length(CLOZUK.duplicated.removed) >= 1){
  CLOZUK_duplicate_SNPS <- CLOZUK.data$SNP[(duplicated(CLOZUK.data$SNP) | duplicated(CLOZUK.data$SNP, fromLast = TRUE)) ]
  CLOZUK.data <- CLOZUK.data[!(duplicated(CLOZUK.data$SNP) | duplicated(CLOZUK.data$SNP, fromLast = TRUE)) ]
}
length(CLOZUK.data$SNP)

cat("Number of SNPS in ",Validation_name," after duplication check Chr:",chromosome.number, "N=" ,nrow(PGC.data.frame))

### Replacing PGC values###
PGC.integers.to.change <- PGC.data.frame[,.I[grep(paste0("^",chromosome.number,":\\d+:\\w+:\\w+"),SNP, perl = T, invert = T)]]
PGC.divisive.names <- PGC.data.frame[PGC.integers.to.change]
PGC.alternative <- PGC.data.frame[,SNP := paste0(CHR,":",BP)]
PGC.data.frame <- PGC.data.frame[,SNP := paste0(CHR,":",BP)]
#PGC.alternative <- PGC.alternative[,c("CHR","SNP","BP","A1","A2","OR"),with = F]
PGC.alternative <- PGC.alternative[, Column_headers, with = F]

if (BETA == FALSE & OR == TRUE){  

### Adding BETA column to data
log.to.odds(PGC.alternative)
PGC.alternative[,BETA := e$PGC.BETA]
PGC.data.frame[,BETA := e$PGC.BETA]
}else{
  stop("neither BETAs or Odds ratios supplied in the Training set datatable, please calculate before merging with Genotype")
}

### Changes the CLOZUK identifiers ###
CLOZUK.integers.to.change <- CLOZUK.data[,.I[grep(paste0("^",chromosome.number,":\\d+:\\w+:\\w+"),SNP, perl = T,invert = T)]]
CLOZUK.divisive.names <- CLOZUK.data[CLOZUK.integers.to.change]
CLOZUK.alternative <- CLOZUK.data[,SNP := paste0(CHR,":",BP)]
CLOZUK.alternative <- CLOZUK.alternative[,c("CHR","SNP","BP","A1","A2"),with = F]

# Simple tracking of each row for future checks
CLOZUK.alternative$place <- c(1:length(CLOZUK.alternative$CHR))
PGC.alternative$place <- c(1:length(PGC.alternative$CHR))

### Merging based on BP position ###
combined.CLOZUK.PGC <- merge(x = PGC.alternative, y = CLOZUK.alternative, by = c('BP','CHR'), all = F, sort = F)

### Remove Duplicated SNPs here ###
combined_a <- which(duplicated(combined.CLOZUK.PGC$SNP.x))
combined_a_rev <- which(duplicated(combined.CLOZUK.PGC$SNP.x, fromLast = T))
combined_b <- which(duplicated(combined.CLOZUK.PGC$SNP.y))
combined_b_rev <- which(duplicated(combined.CLOZUK.PGC$SNP.y,fromLast = T))

### One-line duplication check ###
length(combined.CLOZUK.PGC$SNP.x)
if (all(combined_a == combined_b)){
  duplicate_SNPS <- combined.CLOZUK.PGC$SNP.x[(duplicated(combined.CLOZUK.PGC$SNP.x) | duplicated(combined.CLOZUK.PGC$SNP.x, fromLast = TRUE)) ]
  combined.CLOZUK.PGC <- combined.CLOZUK.PGC[!(duplicated(combined.CLOZUK.PGC$SNP.x) | duplicated(combined.CLOZUK.PGC$SNP.x, fromLast = TRUE)) ]
}else{
  stop("duplicated SNPs in training datset, refer back and correct")
}
length(combined.CLOZUK.PGC$SNP.x)

### Checking and limiting to SNPs with only one allele for each A1 and A2
Checking_length_of_alleles(combined.CLOZUK.PGC)

# Change all alleles to uppercase: saving time method
# below it states:
## assign columns you want to change to object cols
## take all rows of data.table; apply the toupper function to Subset of Data columns named in cols object
cols <- c('A1.x','A2.x','A1.y','A2.y')
combined.CLOZUK.PGC[, (cols) := lapply(.SD, toupper), .SDcols = cols ]

cat("Number of SNPs BEFORE flipping between CLOZUK and PGC Chr:", chromosome.number , "N=", nrow(combined.CLOZUK.PGC))

if (nrow(combined.CLOZUK.PGC) > nrow(PGC.alternative) | nrow(CLOZUK.alternative)){
  warning("combined dataset size is larger than one/both of the input datasets, check for duplicates")
}

a <- which(combined.CLOZUK.PGC$A1.x == "C" & combined.CLOZUK.PGC$A2.x == "G"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.x == "G" & combined.CLOZUK.PGC$A2.x == "C"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.x == "A" & combined.CLOZUK.PGC$A2.x == "T"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.x == "T" & combined.CLOZUK.PGC$A2.x == "A"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}

cat(Validation_name,"_", Training_name, "Chr:",chromosome.number,"remove A-T, C-G Step one: N=" ,nrow(combined.CLOZUK.PGC))

a <- which(combined.CLOZUK.PGC$A1.y=="C" & combined.CLOZUK.PGC$A2.y=="G"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.y=="G" & combined.CLOZUK.PGC$A2.y=="C"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.y=="A" & combined.CLOZUK.PGC$A2.y=="T"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}
a <- which(combined.CLOZUK.PGC$A1.y=="T" & combined.CLOZUK.PGC$A2.y=="A"); if (length(a)>0) {combined.CLOZUK.PGC<-combined.CLOZUK.PGC[-a]}


cat("CLOZUK_PGC Chr:", chromosome.number, "remove A-T, C-G Step two: N=", nrow(combined.CLOZUK.PGC))

## FLIPPING ##
# A1.x/A2.x are Training set alleles
# A1.y/A2.y are Genotype set alleles

a <- which (combined.CLOZUK.PGC$A1.y == combined.CLOZUK.PGC$A2.x & combined.CLOZUK.PGC$A2.y == combined.CLOZUK.PGC$A1.x)
combined.CLOZUK.PGC$BETA[a] <- (-combined.CLOZUK.PGC$BETA[a])

# Change the values of the OR to suit Allele changes #
if (length(a) >= 1 & OR == TRUE) {
  PGC.OR <- combined.CLOZUK.PGC$OR[a]
  change.odds(PGC.OR)
  combined.CLOZUK.PGC$OR[a]<- e$PGC.NEW.OR
  rm(PGC.NEW.OR, envir = e)
  rm(PGC.OR)
}

if (length(a) >= 1){
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
if (length(a) >= 1 & OR == TRUE){
  combined.CLOZUK.PGC$BETA[a] <- (-combined.CLOZUK.PGC$BETA[a])
  PGC.OR <- combined.CLOZUK.PGC$OR[a]
  change.odds(PGC.OR)
  combined.CLOZUK.PGC$OR[a]<- e$PGC.NEW.OR
  rm(PGC.NEW.OR, envir = e)
  rm(PGC.OR)
}
# If OR is not there #
if (length(a) >= 1){
  combined.CLOZUK.PGC$A1.x[a] <- as.character(combined.CLOZUK.PGC$A1.y[a])
  combined.CLOZUK.PGC$A2.x[a] <- as.character(combined.CLOZUK.PGC$A2.y[a])
}

a <- which(combined.CLOZUK.PGC$A1.y != combined.CLOZUK.PGC$A1.x | combined.CLOZUK.PGC$A2.y!=combined.CLOZUK.PGC$A2.x)
if (length(a)>0) combined.CLOZUK.PGC <- combined.CLOZUK.PGC[-a,]
cat("Chr:", chromosome.number,", remove A-T, C-G and other mismatches: N=", nrow(combined.CLOZUK.PGC))

### Find  out where in the original table the SNP identifiers are
PGC.final.index <- combined.CLOZUK.PGC$place.x
PGC.A1.allele.changes <- combined.CLOZUK.PGC$A1.x
PGC.A2.allele.changes <- combined.CLOZUK.PGC$A2.x
if(OR == TRUE){
PGC.OR.changes <- combined.CLOZUK.PGC$OR
}
PGC.BETA.changes <- combined.CLOZUK.PGC$BETA
CLOZUK.final.index <- combined.CLOZUK.PGC$place.x


### Alter original files to have the right allele switching and the right BETA coefficient #
PGC.data.frame <- PGC.data.frame[PGC.final.index,A1:= PGC.A1.allele.changes]
PGC.data.frame <- PGC.data.frame[PGC.final.index,A2:= PGC.A2.allele.changes]
PGC.data.frame <- PGC.data.frame[PGC.final.index,OR:= PGC.OR.changes]
PGC.data.frame <- PGC.data.frame[PGC.final.index,BETA:= PGC.BETA.changes]
oldnames_training <- colnames(PGC.data.frame)

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
if (BETA == FALSE){
namesPGC <- names(PGC.data.frame)
OR_position <- which("OR" == oldnames_training)
BETA_position <- which("BETA" == namesPGC)
namesPGC <- namesPGC[c(1:(OR_position-1),BETA_position,OR_position:length(names_PGC))]
PGC.data.frame <- PGC.data.frame[, namesPGC, with = F]
PGC.data.frame <- PGC.data.frame[, OR := NULL]
}

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
  filename.CLOZUK.together <- paste0("./output/",Validation_name,"_chr", chromosome.number,"_chr.pos.txt")
  new.PGC.table <- paste0("./output/",Training_name,"_table", chromosome.number,"_new.txt")
  filename.common.snps <- paste0("./output/chr", chromosome.number, Training_name,"_", Validation_name,"_common_SNPs.txt")
  filename.duplicate.snps <- paste0("./output/extracted_Duplicate_snps_",Validation_name,"_", Training_name,"_chr",chromosome.number,".txt")
} else {  
  filename.CLOZUK.together <- paste0("./output/", Validation_name,"_chr", chromosome.number,"_chr.pos.txt")
  new.PGC.table <- paste0("./output/", Training_name,"_table", chromosome.number,"_new.txt")
  filename.common.snps <- paste0("./output/chr", chromosome.number, Training_name,"_", Validation_name,"_common_SNPs.txt")
  filename.duplicate.snps <- paste0("./output/extracted_Duplicate_snps_",Validation_name,"_", Training_name,"_chr",chromosome.number,".txt")
}

# Write update file for plink
write.table(CLOZUK_together, file = filename.CLOZUK.together, quote = F, row.names = F, col.names = F)

# Write out new PGC table
write.table(PGC.data.frame, file = new.PGC.table, row.names = F, quote = F)

## Write out the common SNPs between CLOZUK and PGC
write(combined.CLOZUK.PGC$SNP.x, file = filename.common.snps)

### write SNPs back to file
write(duplicate_SNPS, file = filename.duplicate.snps)

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


