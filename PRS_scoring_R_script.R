##################################################
# Checking location for serial or batch analysis #
##################################################

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

### Library
library(data.table)

### environment for functions
e <- new.env()

### Function which assigns genes to the SNP data ####
Assigning_genes <- function(gene.regions, BP.clumped.SNPs, clumped_SNPs, outputfilename){
GR <- deparse(substitute(gene.regions))
  
if(GR == "MAGMA.gene.regions.for.chromosome"){
  
  test_vector <- rep(0,nrow(clumped_SNPs))
  gene_ID_vector <- rep(0,nrow(clumped_SNPs))
  
    for (i in 1:nrow(clumped_SNPs)) {
    tmp_gene <- which ((gene.regions$BP_1 <= BP.clumped.SNPs[i] & gene.regions$BP_2 >= BP.clumped.SNPs[i]) == T)
    if (length(tmp_gene) == 1) {
      test_vector[i] <- gene.regions$Gene_symbol[tmp_gene]
      gene_ID_vector[i] <- gene.regions$ID[tmp_gene]
    }else{
      test_vector[i] <- NA
      gene_ID_vector[i] <- NA
    }
    }
  ## This is specialised notation for a data.table R object
  ## it is essentially saying, take all rows (defalt blank before comma), assign column Genes to the R object test_vector
  ## data.table will automatically create a new column because the column Gene does not exist in the original data.table
  
  clumped_SNPs[, Gene_name := test_vector]
  clumped_SNPs[, Gene_ID := gene_ID_vector]
  
  ## Which SNPs are not in genes?
  ## .I prints a vector of the rows that match the expression within the square brackets
  ## .() would do a similar thing but keeps the format of the object as a data.table
  index.non.genic.SNPs <- clumped_SNPs[,.I[which(is.na(Gene_name))]]
  
  ## Remove those SNPs from the table
  Gene_clumped_SNPs <- clumped_SNPs[-index.non.genic.SNPs]
  Gene_clumped_SNPs <- Gene_clumped_SNPs[order(BP, decreasing = F)]
  
  ## Print out new table (replace with fwrite once you update it)
  write.table(Gene_clumped_SNPs, 
              file = paste0("./output/", outputfilename, chromosome.number, ".txt"), 
              quote = F, 
              row.names = F)
  assign("Gene_clumped_SNPs",Gene_clumped_SNPs, envir = e)
  }

if(GR == "MAGMA.expanded.gene.regions.for.chromosome"){
  
  test_vector <- rep(0,nrow(clumped_SNPs))
  gene_ID <- rep(0,nrow(clumped_SNPs))
  
  for (i in 1:nrow(clumped_SNPs)) {
    tmp_gene <- which ((gene.regions$BP_1 <= BP.clumped.SNPs[i] & gene.regions$BP_2 >= BP.clumped.SNPs[i]) == T)
    if (length(tmp_gene) == 1) {
      test_vector[i] <- gene.regions$Gene_symbol[tmp_gene]
      gene_ID[i] <- gene.regions$ID[tmp_gene]
    }else{
      test_vector[i] <- NA
      gene_ID[i] <- NA
    }
  }
  
## This is specialised notation for a data.table R object
## it is essentially saying, take all rows (defalt blank before comma), assign column Genes to the R object test_vector
## data.table will automatically create a new column because the column Gene does not exist in the original data.table

  clumped_SNPs[, Gene_name := test_vector]
  clumped_SNPs[, Gene_ID := gene_ID]

## Which SNPs are not in genes?
## .I prints a vector of the rows that match the expression within the square brackets
## .() would do a similar thing but keeps the format of the object as a data.table
index.non.genic.SNPs <- clumped_SNPs[,.I[which(is.na(Gene_name))]]

## Remove those SNPs from the table
Gene_clumped_SNPs <- clumped_SNPs[-index.non.genic.SNPs]
Gene_clumped_SNPs <- Gene_clumped_SNPs[order(BP, decreasing = F)]

## Print out new table (replace with fwrite once you update it)
write.table(Gene_clumped_SNPs, 
            file = paste0("./output/", outputfilename, chromosome.number, ".txt"), 
            quote = F, 
            row.names = F)

assign("Expanded_Gene_clumped_SNPs",Gene_clumped_SNPs, envir = e)
}
}


## PERFORM CLUMPING (not added yet) ##
## Read in clumped data
clumped_SNPs <- fread(paste0("./output/CLOZUK_PGC_CLUMPED_FINAL_chr",chromosome.number,".txt"))
wd <-getwd()

## read in MAGMA's the gene regions 
MAGMA.gene.regions <- fread("NCBI37.3.gene.loc",colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
colnames(MAGMA.gene.regions) <- c("ID","Chromosome","BP_1","BP_2","strand","Gene_symbol")

## Limit the data down to the specific chromsome the array is on
Index.for.chromosome <- MAGMA.gene.regions[,.I[grep(paste0("^",chromosome.number),Chromosome, perl = T, invert = F)]]
MAGMA.gene.regions.for.chromosome <- MAGMA.gene.regions[Index.for.chromosome]
MAGMA.gene.regions.for.chromosome <- MAGMA.gene.regions.for.chromosome[,c('ID','Chromosome','BP_1','BP_2','Gene_symbol'),with = F]
########## DO I NEED THE BELOW?#########
write.table(MAGMA.gene.regions.for.chromosome,file = paste0("Chromosome_", chromosome.number, "_MAGMA_GENE_REGIONS_NCBI37.3.txt"),quote = F,row.names = F,col.names = F)
########################################
## read in -35kb upstream and 10kb downstream regions for genes

# Note that there is an overlap between genes
# This shouldn't make any difference at the moment but if you were to specify which genes are used as the p-values there might be a problem determining in which areas the SNP's belong to
# Additionally, you should find out how magma does it's analysis

##### different script: keep separate #####
MAGMA.expanded.gene.regions <- MAGMA_different_gene_regions <- fread("NCBI37_-35kb_+10kb.gene.loc",colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
colnames(MAGMA.expanded.gene.regions) <- c("ID","Chromosome","BP_1","BP_2","strand","Gene_symbol")

## Limit the data down to the specific chromsome the array is on
Index.for.chromosome.expanded <- MAGMA.expanded.gene.regions[,.I[grep(paste0("^",chromosome.number),Chromosome, perl = T, invert = F)]]
MAGMA.expanded.gene.regions.for.chromosome <- MAGMA.expanded.gene.regions[Index.for.chromosome]
MAGMA.expanded.gene.regions.for.chromosome <- MAGMA.expanded.gene.regions.for.chromosome[,c('ID','Chromosome','BP_1','BP_2','Gene_symbol'), with = F]
write.table(MAGMA.expanded.gene.regions.for.chromosome, file = paste0("Chromosome_", chromosome.number, "_MAGMA_GENE_REGIONS_EXPANDED_NCBI37.3.txt"), quote = F, row.names = F, col.names = F)
###################################
## remove orginal tables
rm(MAGMA.gene.regions)

###################################
rm(MAGMA.expanded.gene.regions)
###################################

## Using Data tables to get the gene name printed in the row against the SNP
## Mapping Genes to SNP's and storing in Data.table

BP.clumped.SNPs <- clumped_SNPs$BP

Assigning_genes(MAGMA.gene.regions.for.chromosome, BP.clumped.SNPs = BP.clumped.SNPs, clumped_SNPs = clumped_SNPs, outputfilename = "Genomic_1000kb_r2_zeropoint2_PGC_CLOZUK_chr_")
Assigning_genes(MAGMA.expanded.gene.regions.for.chromosome, BP.clumped.SNPs = BP.clumped.SNPs, clumped_SNPs = clumped_SNPs, outputfilename = "Genomic_expanded_1000kb_r2_zeropoint2_PGC_CLOZUK_chr_")

GENES_to_snps <- scan(file = paste0("./output/CLOZUK_PRS_CLUMPED_chr",chromosome.number,".genes.annot"), what = "", sep = "\n")
y <- strsplit(GENES_to_snps, "[[:space:]]+")
names(y) <- sapply(y, '[[', 1)
y <- lapply(y, '[', -1)

y[[1]] <- NULL
y[[1]] <- NULL

adding_unread_genes <- function(Gene_clumped_SNPs, MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number){
for (i in 1:length(y)) {
  index_of_gene <- which(Gene_clumped_SNPs$Gene_ID == names(y[i]))
  if (length(index_of_gene) == 0){
    gene.to.include <- names(y[i])
    SNPs.to.include <- y[[i]][2:length(y[[i]])]
    
    Gene.names.for.table <- which(gene.to.include == MAGMA.gene.regions.for.chromosome$ID)
      for (f in 1:length(SNPs.to.include)) {
        Base.pairs.index <- which(clumped_SNPs$SNP == SNPs.to.include[f])
        new.row <- list(chromosome.number, 
                        SNPs.to.include[f], 
                        clumped_SNPs$BP[Base.pairs.index], 
                        clumped_SNPs$P[Base.pairs.index], 
                        MAGMA.gene.regions.for.chromosome$Gene_symbol[Gene.names.for.table], 
                        gene.to.include)
        Gene_clumped_SNPs <- rbindlist(list(Gene_clumped_SNPs,new.row))
      }
  }else{
  gene.names.in.y <- Gene_clumped_SNPs[index_of_gene]
  unrecorded_SNPs <- which(!y[[i]][2:length(y[[i]])] %in% gene.names.in.y$SNP)
  
  if (length(unrecorded_SNPs) != 0 ) {
    
    gene.to.include <- names(y[i])
    Gene.names.for.table <- which(gene.to.include == MAGMA.gene.regions.for.chromosome$ID)
    SNPs.to.include <- y[[i]][2:length(y[[i]])]
    
    for (f in 1:length(unrecorded_SNPs)) {
    Base.pairs.index <- which(clumped_SNPs$SNP == SNPs.to.include[unrecorded_SNPs[f]])
    new.row <- list(chromosome.number, 
                    SNPs.to.include[unrecorded_SNPs[f]], 
                    clumped_SNPs$BP[Base.pairs.index], 
                    clumped_SNPs$P[Base.pairs.index], 
                    MAGMA.gene.regions.for.chromosome$Gene_symbol[Gene.names.for.table], 
                    MAGMA.gene.regions.for.chromosome$Gene_symbol[Gene.names.for.table],
                    gene.to.include)
    Gene_clumped_SNPs <- rbindlist(list(Gene_clumped_SNPs,new.row))
  }
  } #name of gene not in list then add row of SNPs with gene identifier to the table()
}
}
  assign("test_data_frame",Gene_clumped_SNPs,envir = .GlobalEnv)
}

adding_unread_genes(e$test_clumped_snps,MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number)  
names(test_data_frame) <- c("CHR" ,"SNP", "BP", "P", "Gene_name", "Gene_ID","BP_1","BP_2")

#ALTER HERE TO MATCH PYTHON# 
#ONE for the gene names and one for the gene number

#!!!!!!!!!!!!!!!! Check if the SNPs are being assigned to the correct genes, if not then you will ignore some genes in the analysis  
which(duplicated(test_data_frame$SNP,fromLast = T))

write.table(test_data_frame,file = paste0("./output/MAGMA_Gene_regions_for_python_script_chr_",chromosome.number,".txt")
##############

##Read in PGC.data

PGC <- fread("PGC_table22.txt")

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
  
p.value.thresholds <- c(0.0001,0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5)

scoring_output_file <- paste0("scoring_PGC_CLOZUK_chromosome_", chromosome.number)

## check that it is merging properly here (after analysis is run)
combined.CLOZUK.PGC.clumped.Genomic.SNPs <- merge(test_data_frame, PGC, by.x="SNP", by.y="SNP", all=F, sort=F)
combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1 <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1)
combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2 <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2)

Genes_index <- unique(combined.CLOZUK.PGC.clumped.Genomic.SNPs$Gene_ID)
Genes_used <- matrix(ncol = 2)

for (i in 1:length(p.value.thresholds)) {
  for (l in 1:length(Genes_index)) {

  a <- copy(combined.CLOZUK.PGC.clumped.Genomic.SNPs)    
  SNPs <- a[, .I[which(Gene_ID == Genes_index[l] & P.x <= p.value.thresholds[i])]]
  
    if (length(SNPs) != 0){
    a <- a[SNPs, .(SNP,A1,BETA)]
  
    filename <- paste0('score_1/chr22_test_', Genes_index[l],'_', p.value.thresholds[i],".score")
  

    write.table(file = filename, a, row.names = F, col.names = F, quote = F, sep="\t")
  
    command <- paste0('plink --bfile CLOZUK_GWAS_BGE_CLUMPED_chr22 --score ', filename, " --out /Users/JJ/Documents/CLOZUK_PRS/Profiles/chr22_1/chr22_test_", Genes_index[l], "_", p.value.thresholds[i], "_a")
    system(command)
    rm(a)
    genes_and_pvalue <- c(Genes_index[l],p.value.thresholds[i])
    Genes_used <- rbind(Genes_used,genes_and_pvalue)
    }else{
      next()
    }
  }
}
write.table(Genes_used, file = "Index_of_genes_and_pval_1.txt", quote = F, col.names = T, row.names = F)
  #score
  
for (i in 1:(length(p.value.thresholds)))
{
# format the p.values so that exponential notation is NOT used, could result in confusion in BASH scripting
  scoring.output.filename <- paste0(scoring_output_file, "_", format(p.value.thresholds[i], scientific =  F), ".score")
  SNPs.lower.than.threshold <- which(combined.CLOZUK.PGC.clumped.Genomic.SNPs$P.x <= p.value.thresholds[i]) 
  
  if (length(SNPs.lower.than.threshold ) > 0) {
    tmp <- combined.CLOZUK.PGC.clumped.Genomic.SNPs[SNPs.lower.than.threshold, ]
    final.output <- tmp[, c("SNP","A1","BETA"), with = F]
  }
  
  if (nrow(out)>0) {
    write.table(file = scoring.output.filename, final.output, row.names = F, col.names = F, quote = F, sep="\t")
  }
}
## END script!
#######################################
