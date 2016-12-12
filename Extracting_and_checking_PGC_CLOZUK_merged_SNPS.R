### Finding list of useful SNPs from both CLOZUK and PGC ### 

### Chromosome 22 ###

# Add new environment for variables created within a function
e <- new.env()

# ### Checking to see if two columns have alleles swapped ###
Checking_allele_swapping <- function(alteredtable1,table2,which.is.combined = c("NONE","PGC","CLOZUK")){
  if (which.is.combined == "PGC") {
      Allele.x <- alteredtable1[,c("A1.x","A2.x"),with = F]
      Allele.y <- table2[,c("A1","A2"),with = F]
  }
  if(which.is.combined == "CLOZUK") {
    Allele.x <- alteredtable1[,c("A1.y","A2.y"),with = F]
    Allele.y <- table2[,c("A1","A2"),with = F]
  }
  if(which.is.combined == "NONE"){
  Allele.x <- alteredtable1[,c("A1","A2"),with = F]
  Allele.y <- table2[,c("A1","A2"),with = F]
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
########




# Preparing to run in local
chromosome.number <- 22
system_information<-Sys.info()
if (system_information[1] == "Windows") fpath <-  "/Users/JJ/" else fpath <-"/Users/johnhubert/"

# load libraries
library(data.table)

# Set Working directory
setwd(paste0(fpath,"Dropbox/PGC_CLOZUK_GWAS_INPUT/"))

### Adding in PGC data ###
PGC.test.data.frame <- fread(paste0("gzip -dc PGC_table",chromosome.number,".txt.gz"))
cat("Chr:",chromosome.number, nrow(PGC.test.data.frame))

### Read in the CLOZUK data ###
untar(paste0("CLOZUK_GWAS_BGE_chr",chromosome.number,".tar.gz"),files = paste0("CLOZUK_GWAS_BGE_chr",chromosome.number,".bim"))
CLOZUK.data <- fread(paste0("CLOZUK_GWAS_BGE_chr",chromosome.number,".bim"))
cat("Chr:",chromosome.number, nrow(CLOZUK.data))

### Replace the column names ###
setnames(CLOZUK.data, c("CHR","SNP","GENEDIST","BP","A1","A2"))

# Set Working directory
setwd(paste0(fpath,"Documents/PR54/PGC_CLOZUK_PRS/"))

### Adding in combined data ###
Combined.PGC.CLOZUK <- fread(paste0("gzip -dc output/PGC_CLOZUK_SNP_table",chromosome.number,".txt.gz"))
cat("Chr:",chromosome.number, nrow(Combined.PGC.CLOZUK))

### Adding in extra info ###
CLOZUK.altered <- fread(paste0("gzip -dc extrainfo/CLOZUK_altered_names_chr",chromosome.number,".txt.gz"))
PGC.altered <- fread(paste0("gzip -dc extrainfo/PGC_altered_names_chr",chromosome.number,".txt.gz"))
PGC.index <- scan(paste0("extrainfo/PGC.merged.index.chr",chromosome.number,".txt"), sep = " ")
CLOZUK.index <- scan(paste0("extrainfo/CLOZUK.merged.index.chr",chromosome.number,".txt"), sep = " ")

### checking everything is above board ###
Checking_allele_swapping(PGC.altered,PGC.test.data.frame,which.is.combined = "NONE")
Checking_allele_swapping(CLOZUK.altered,CLOZUK.data,which.is.combined = "NONE")

if (length(e$CLOZUK.alteredflipped.alleles) != 0) {
  stop("alleles have been flipped on CLOZUK")
}

checking.which.OR.do.not.equal.each.other <- which(PGC.altered$OR != PGC.test.data.frame$OR)
problematic.SNPs <- match(e$PGC.alteredflipped.alleles, checking.which.OR.do.not.equal.each.other)  
problematic.SNPs <- which(is.na(problematic.SNPs))
check1_index <- e$PGC.alteredflipped.alleles[problematic.SNPs]
check1 <- PGC.altered[check1_index]
check2 <- PGC.test.data.frame[check1_index]

if (all(check1$OR == 1) & all(check2$OR == 1)){
  cat("You have",length(check1$CHR),"flipped SNPs with OR = 1")
}else{
  stop("There is an uneven amount of flipped SNPs")
}

### Checking conversions of SNPs ###
checking.duplications.PGC <- which (duplicated(Combined.PGC.CLOZUK$SNP.x))
checking.duplications.CZK <- which (duplicated(Combined.PGC.CLOZUK$SNP.y))

if (length(checking.duplications.CZK) != 0 & length(checking.duplications.PGC) != 0) {
  stop("There are duplicated SNPs common between CLOZUK and PGC")
}

### Extracting SNPs for analysis ###


### Taring multiple files together ###

### TODO ###

# Incorporate the new ID's into CLOZUK data including the genotype information etc... tar all together and keep in bim format#
# Extract a list of common SNP's for each chromosome as a vector and print out to file.txt
# print out new PGC table that has altered OR and BETA to be used for analysis that is usable



  