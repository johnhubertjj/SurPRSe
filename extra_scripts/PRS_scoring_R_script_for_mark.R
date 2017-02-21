# NOTE for Mark
### currently runs for Chromosome 22 only,and only runs locally at the moment, aiming to set up for job array (using array variable as chromosome number)###

### Library

library(data.table)

### Read in updated_raw data
setwd("/Users/johnhubert/Dropbox/PhD_clumping/")
PGC <- fread("PGC_table22.txt")
CLOZUK <- fread ("CLOZUK_GWAS_BGE_chr22.bim")
Common_SNPs <- scan("chr22PGC_CLOZUK_common_SNPs.txt",what = "numeric")

### Check for duplicated SNPs.PGC
PGCa <- which(duplicated(PGC$SNP))
PGCb <- which(duplicated(PGC$SNP,fromLast = T))
PGC_duplicated_SNPs <- c(a,b)
removed_SNPs <- PGC[PGC_duplicated_SNPs]
removed_SNPs <- removed_SNPs$SNP

### Check for duplicated SNPs.CLOZUK
### Add rownames so it is more readable

cloza <- which(duplicated(CLOZUK$V2))
clozb <- which(duplicated(PGC$SNP,fromLast = T))
removed_SNPs_cloz <- c(cloza,clozb)
removed_SNPs_cloz <- CLOZUK[removed_SNPs_cloz]
removed_SNPs_cloz <- removed_SNPs_cloz$V2

### Extract duplicate SNPs
extract_SNPs <-c(removed_SNPs,removed_SNPs_cloz)

### write SNPs back to file
write(extract_SNPs,file = "extracted_Duplicate_snps.txt")
which(removed_SNPs %in% Common_SNPs)


# Add a system command that runs plink and add some dependencies so that R will run #

# Checking that the relevant file exists

while (!file.exists(your.file.name)) {
  Sys.sleep(1)
}



## PERFORM CLUMPING (not added yet) ##
## Change from output being tmp1.txt ##

## Read in clumped data
clumped_SNPs <- fread("tmp1.txt")
wd <-getwd()
chromosome.number <- 22

## read in MAGMA's the gene regions 
MAGMA.gene.regions <- fread("NCBI37.3.gene.loc",colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
Example.gene.regions <- fread("example_clumping_gene_regions.txt")

MAGMA.gene.regions22 <- MAGMA.gene.regions[,.I[grep(paste0("^",chromosome.number),V2, perl = T,invert = F)]]
MAGMA.gene.regions22 <- MAGMA.gene.regions[MAGMA.gene.regions22]
MAGMA.gene.regions22 <- MAGMA.gene.regions22[,c('V2','V3','V4','V6'),with = F]
write.table(MAGMA.gene.regions22,file = "Chromosome_22_MAGMA_GENE_REGIONS_NCBI37.3.txt",quote = F,row.names = F,col.names = F)



## Using Data tables to get the gene name printed in the row against the SNP
## Mapping Genes to SNP's and storing in Data.table

BP.clumped.SNPs <- clumped_SNPs$BP
gene.location <- rep(0,length(BP.clumped.SNPs))
for (i in 1:nrow(clumped_SNPs)) {
  tmp_gene <- which ((MAGMA.gene.regions22$V3 <= BP.clumped.SNPs[i] & MAGMA.gene.regions22$V4 >= BP.clumped.SNPs[i]) == T)
  if (length(tmp_gene) == 1) {
    test_vector[i] <- MAGMA.gene.regions22$V6[tmp_gene]
  }else{
    test_vector[i] <- NA
  }
}

## This is specialised notation for a data.table R object
## it is essentially saying, take all rows (defalt blank before comma), assign column Genes to the R object test_vector
## data.table will automatically create a new column because the column Gene does not exist in the original data.table

clumped_SNPs[, Gene := test_vector]


## Which SNPs are not in genes?
## .I prints a vector of the rows that match the expression within the square brackets
## .() would do a similar thing but keeps the format of the object as a data.table
index.non.genic.SNPs <- clumped_SNPs[,.I[which(is.na(Gene))]]

## Remove those SNPs from the table
Gene_clumped_SNPs <- clumped_SNPs[-index.non.genic.SNPs]
Gene_clumped_SNPs <- Gene_clumped_SNPs[order(BP, decreasing = F)]

## Print out new table (replace with fwrite once you update it)
write.table(Gene_clumped_SNPs, file = "/Users/johnhubert/Desktop/Genomic_1000kb_r2_zeropoint2_PGC_CLOZUK.txt",quote = F,row.names = F)



## Plink --gene-report tmp1.txt Chromosome_22_MAGMA_GENE_REGIONS_NCBI37.3.txt
##       --out CLOZUK_BGE_GENE_REGIONS

# Plink --bfile CLOZUK_GWAS_BGE_chr22
#       --make-set Chromosome_22_MAGMA_GENE_REGIONS_NCBI37.3.txt
#       --out testlol
#       --write-set



#Scoring per gene: rough and undescriptive at the moment

# SCORING:
  
sig<-c(1e-4,1e-3,1e-2,0.05,0.1,0.2,0.3,0.4,0.5)

scoring_output_file <-"scoring_PGC_CLOZUK_chromosome_22.txt"

## check that it is merging properly here
combined.CLOZUK.PGC.clumped.Genomic.SNPs <- merge(Gene_clumped_SNPs, PGC, by.x="SNP", by.y="SNP", all=F, sort=F)
m$A1 <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1)
m$A2 <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2)

#score
for (i in 1:(length(sig)))
{
  fn <- sprintf("%s_%f.score", scoring_output_file, sig[i])
  a <- which(m$P.x<=sig[i]) 
  if (length(a)>0)
  {
    tmp<-m[a,]
    out<-tmp[,c("SNP","A1","BETA"), with = F]
  }
  if (nrow(out)>0) {write.table(file=fn, out, row.names=F, col.names=F, quote=F, sep="\t")}
}

# PROFILING

sig<-c(1e-4,1e-3,1e-2,0.05,0.1,0.2,0.3,0.4,0.5)

for (i in 1:(length(sig)))
{
  for (chr in 22)
  {
    command_txt<-sprintf("plink2 --silent --bfile  CLOZUK_GWAS_BGE_chr22 --score %s_%f.score,  --out  profiles/chr%i.clump_r0.2_1000kb_%f", chr, sig[i], chr, sig[i])
    
    system('bash -l',input=c("shopt -s expand_aliases",command_txt))
  }
} #sig


