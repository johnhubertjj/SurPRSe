### Finding list of useful SNPs from both CLOZUK and PGC ### 

### Chromosome 22 ###

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

# Set Working directory
setwd(paste0(fpath,"Documents/PR54/PGC_CLOZUK_PRS/"))

### Adding in combined data ###
Combined.PGC.CLOZUK <- fread(paste0("gzip -dc output/PGC_CLOZUK_SNP_table",chromosome.number,".txt.gz"))
cat("Chr:",chromosome.number, nrow(Combined.PGC.CLOZUK))

### Adding in extra info ###
CLOZUK.altered <- fread(paste0("gzip -dc extrainfo/CLOZUK_altered_names_chr",chromosome.number,".txt.gz"))
PGC.altered <- CLOZUK.altered <- fread(paste0("gzip -dc extrainfo/PGC_altered_names_chr",chromosome.number,".txt.gz"))
PGC.index <- scan(paste0("extrainfo/PGC.merged.index.chr",chromosome.number,".txt"), sep = " ")
CLOZUK.index <- scan(paste0("extrainfo/CLOZUK.merged.index.chr",chromosome.number,".txt"), sep = " ")



