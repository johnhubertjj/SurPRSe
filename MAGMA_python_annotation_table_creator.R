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
    BP1 <- rep(0,nrow(clumped_SNPs))
    BP2 <- rep(0,nrow(clumped_SNPs))
    
    for (i in 1:nrow(clumped_SNPs)) {
      tmp_gene <- which ((gene.regions$BP_1 <= BP.clumped.SNPs[i] & gene.regions$BP_2 >= BP.clumped.SNPs[i]) == T)
      if (length(tmp_gene) == 1) {
        test_vector[i] <- gene.regions$Gene_symbol[tmp_gene]
        gene_ID_vector[i] <- gene.regions$ID[tmp_gene]
        BP1 <- gene.regions$BP1[tmp_gene]
        BP2 <- gene.regions$BP2[tmp_gene]
        
      }else{
        test_vector[i] <- NA
        gene_ID_vector[i] <- NA
        BP1[i] <- NA
        BP2[i] <- NA
      }
    }
    ## This is specialised notation for a data.table R object
    ## it is essentially saying, take all rows (defalt blank before comma), assign column Genes to the R object test_vector
    ## data.table will automatically create a new column because the column Gene does not exist in the original data.table
    
    clumped_SNPs[, Gene_name := test_vector]
    clumped_SNPs[, Gene_ID := gene_ID_vector]
    clumped_SNPs[, BP_1 := BP1]
    clumped_SNPs[, BP_2 := BP2]
    
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
    BP1 <- rep(0,nrow(clumped_SNPs))
    BP2 <- rep(0,nrow(clumped_SNPs))
    
    for (i in 1:nrow(clumped_SNPs)) {
      
      tmp_gene <- which ((gene.regions$BP_1 <= BP.clumped.SNPs[i] & gene.regions$BP_2 >= BP.clumped.SNPs[i]) == T)
      
      if (length(tmp_gene) == 1) {
        test_vector[i] <- gene.regions$Gene_symbol[tmp_gene]
        gene_ID[i] <- gene.regions$ID[tmp_gene]
        BP1 <- gene.regions$BP1[tmp_gene]
        BP2 <- gene.regions$BP2[tmp_gene]
        
      }else{
        test_vector[i] <- NA
        gene_ID[i] <- NA
        BP1[i] <- NA
        BP2[i] <- NA
      }
    }
    
    ## This is specialised notation for a data.table R object
    ## it is essentially saying, take all rows (defalt blank before comma), assign column Genes to the R object test_vector
    ## data.table will automatically create a new column because the column Gene does not exist in the original data.table
    
    clumped_SNPs[, Gene_name := test_vector]
    clumped_SNPs[, Gene_ID := gene_ID]
    clumped_SNPs[, BP_1 := BP1]
    clumped_SNPs[, BP_2 := BP2]
    
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

### Function which adds duplicate SNPs which happen to be inside other genes
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
                        MAGMA.gene.regions.for.chromosome$BP_1[Gene.names.for.table],
                        MAGMA.gene.regions.for.chromosome$BP_2[Gene.names.for.table],
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
                          MAGMA.gene.regions.for.chromosome$BP_1[Gene.names.for.table],
                          MAGMA.gene.regions.for.chromosome$BP_2[Gene.names.for.table],
                          gene.to.include)
          Gene_clumped_SNPs <- rbindlist(list(Gene_clumped_SNPs,new.row))
        }
      } #name of gene not in list then add row of SNPs with gene identifier to the table()
    }
  }
  assign("test_data_frame",Gene_clumped_SNPs,envir = .GlobalEnv)
}

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

## remove orginal tables
rm(MAGMA.gene.regions)

# Find clumped SNPs
BP.clumped.SNPs <- clumped_SNPs$BP

# Create table and write to file for reference
Assigning_genes(MAGMA.gene.regions.for.chromosome, 
                BP.clumped.SNPs = BP.clumped.SNPs, 
                clumped_SNPs = clumped_SNPs, 
                outputfilename = "Genomic_1000kb_r2_zeropoint2_PGC_CLOZUK_chr_")

# read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
GENES_to_snps <- scan(file = paste0("./output/CLOZUK_PRS_CLUMPED_chr",chromosome.number,".genes.annot"), what = "", sep = "\n")
y <- strsplit(GENES_to_snps, "[[:space:]]+")
names(y) <- sapply(y, '[[', 1)
y <- lapply(y, '[', -1)

y[[1]] <- NULL
y[[1]] <- NULL

adding_unread_genes(e$test_clumped_snps,MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number)  

names(test_data_frame) <- c("CHR" ,"SNP", "BP", "P", "Gene_name","BP_START","BP_END","GENE")
setcolorder(test_data_frame, c("CHR","SNP","BP","GENE","BP_START","BP_END","P","Gene_name"))

which(duplicated(test_data_frame$SNP,fromLast = T))

write.table(test_data_frame[, c(1:6), with = F],file = paste0("./output/MAGMA_Gene_regions_for_python_script_chr_",chromosome.number,".txt"), quote = F, row.names = F)
write(unique(test_data_frame$Gene_name),file = paste0("./output/PGC_CLOZUK_unique_genes_chr_",chromosome.number,".txt"))


