#! /bin/env/R
# arguments = sex chromosomes involved, serial or batch run, pathway or whole genome?

### Start Timer
ptm <- proc.time()

##################################################
# Checking location for serial or batch analysis #
##################################################

System_information <- Sys.info()
whereami <- System_information['user']
args <- commandArgs(trailingOnly = T)

Training_name <- args[3]
Validation_name <- args[4]
UCSC_file <- args[5]
Validation_dataset <- args[6]
NCBI_file <- args[7]
Gene_regions_to_analyse <- args[8]

# Read in chromosome numbers as an array variable from BASH
chromosomes_to_analyse <- as.numeric(args[c(9:length(args))])
print(chromosomes_to_analyse)

#chromosomes_to_analyse <- c(1:22)
#Gene_regions_to_analyse <- "both"

# Set up library
library("data.table")
library("gtools")

# Set working directory
setwd(".")

###############################
## Location Finding Function ##
###############################

# Find out whether you want to run the script in a batch format or in serial

### environment for functions
e <- new.env()

###################################################
## Function which assigns genes to the SNP data ##
###################################################
Assigning_genes <- function(MAGMA_file, BP.clumped.SNPs, clumped_SNPs, outputfilename, gene.regions = c("normal", "extended", "both"), chromosome.number = l){
  if(gene.regions == "normal" | gene.regions == "both"){
    
    test_vector <- rep(0,nrow(clumped_SNPs))
    gene_ID_vector <- rep(0,nrow(clumped_SNPs))
    BP1 <- rep(0,nrow(clumped_SNPs))
    BP2 <- rep(0,nrow(clumped_SNPs))
    
    for (i in 1:nrow(clumped_SNPs)) {
      tmp_gene <- which ((MAGMA_file$BP_START <= BP.clumped.SNPs[i] & MAGMA_file$BP_END >= BP.clumped.SNPs[i]) == T)
      if (length(tmp_gene) == 1) {
        test_vector[i] <- MAGMA_file$GENE_NAME[tmp_gene]
        gene_ID_vector[i] <- MAGMA_file$Gene[tmp_gene]
        BP1[i] <- MAGMA_file$BP_START[tmp_gene]
        BP2[i] <- MAGMA_file$BP_END[tmp_gene]
        
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
                file = paste0("./output/normal_gene_regions", outputfilename,"chromosome_", chromosome.number,".txt"), 
                quote = F, 
                row.names = F)
    assign(paste0("Gene_clumped_SNPs_chromosome",l), Gene_clumped_SNPs, envir = e)
  }
  
  if(gene.regions == "extended" | gene.regions == "both"){
    test_vector <- rep(0,nrow(clumped_SNPs))
    gene_ID <- rep(0,nrow(clumped_SNPs))
    BP1 <- rep(0,nrow(clumped_SNPs))
    BP2 <- rep(0,nrow(clumped_SNPs))
    
    for (i in 1:nrow(clumped_SNPs)) {
      tmp_gene <- which ((MAGMA_file$BP_start_extended <= BP.clumped.SNPs[i] & MAGMA_file$BP_end_extended >= BP.clumped.SNPs[i]) == T)
      
      if (length(tmp_gene) == 1) {
        test_vector[i] <- MAGMA_file$GENE_NAME[tmp_gene]
        gene_ID[i] <- MAGMA_file$Gene[tmp_gene]
        BP1[i] <- MAGMA_file$BP_start_extended[tmp_gene]
        BP2[i] <- MAGMA_file$BP_end_extended[tmp_gene]
        
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
    clumped_SNPs[, GENE_NAME := test_vector]
    clumped_SNPs[, Gene := gene_ID]
    clumped_SNPs[, BP_start_extended := BP1]
    clumped_SNPs[, BP_end_extended := BP2]
    
    ## Which SNPs are not in genes?
    ## .I prints a vector of the rows that match the expression within the square brackets
    ## .() would do a similar thing but keeps the format of the object as a data.table
    index.non.genic.SNPs <- clumped_SNPs[,.I[which(is.na(GENE_NAME))]]
    
    ## Remove those SNPs from the table
    Gene_clumped_SNPs <- clumped_SNPs[-index.non.genic.SNPs]
    Gene_clumped_SNPs <- Gene_clumped_SNPs[order(BP, decreasing = F)]
    
    ## Print out new table (replace with fwrite once you update it)
    write.table(Gene_clumped_SNPs, 
                file = paste0("./output/Extended_gene_regions_", outputfilename,"chromosome_",chromosome.number ,".txt"), 
                quote = F, 
                row.names = F)
    
    assign(paste0("Expanded_Gene_clumped_SNPs_chromosome_",chromosome.number), Gene_clumped_SNPs, envir = e)
  }
}

###############################################################################
### Function which adds duplicate SNPs which happen to be inside other genes ##
###############################################################################
adding_unread_genes <- function(Gene_clumped_SNPs, MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number, length_of_genes){
# y is the input gene list from MAGMA
# Gene_clumped_SNPs is the updated bim table from the initial assignment of SNPs
if (length_of_genes == "normal"){
  
  for (i in 1:length(y)) {
    # check whether the gene exists in the updated Gene_clumped_SNPs file
    # If it does not, this means that there is a gene inside a gene or the SNPs were only found inside one of the genes 
    index_of_gene <- which(Gene_clumped_SNPs$Gene_ID == names(y[i]))
    if (length(index_of_gene) == 0){
      gene.to.include <- names(y[i])
      # the first elemnt of y[[i]] is the bp start and end so need to take it into account when finding the SNP in question
      SNPs.to.include <- y[[i]][2:length(y[[i]])]
      
      Gene.names.for.table <- which(gene.to.include == MAGMA.gene.regions.for.chromosome$Gene)
      for (f in 1:length(SNPs.to.include)) {
        # clumped_SNPS is the original unaltered Bim file for each chromosome and pathway
        Base.pairs.index <- which(clumped_SNPs$SNP == SNPs.to.include[f])
        # add a new row corresponding to the Gene_clumped_SNPs format
        new.row <- list(MAGMA.gene.regions.for.chromosome$CHR[Gene.names.for.table], 
                        SNPs.to.include[f], 
                        clumped_SNPs$GD[Base.pairs.index],
                        clumped_SNPs$BP[Base.pairs.index],
                        clumped_SNPs$A1[Base.pairs.index],
                        clumped_SNPs$A2[Base.pairs.index],
                        MAGMA.gene.regions.for.chromosome$GENE_NAME[Gene.names.for.table],
                        gene.to.include,
                        MAGMA.gene.regions.for.chromosome$BP_START[Gene.names.for.table],
                        MAGMA.gene.regions.for.chromosome$BP_END[Gene.names.for.table])
        Gene_clumped_SNPs <- rbindlist(list(Gene_clumped_SNPs,new.row))
      }
    }else{
      # check if some SNPs are missing
      gene.names.in.y <- Gene_clumped_SNPs[index_of_gene]
      unrecorded_SNPs <- which(!y[[i]][2:length(y[[i]])] %in% gene.names.in.y$SNP)
      # if they are, identify them and do the same for if the gene was missing
      if (length(unrecorded_SNPs) != 0 ) {
        gene.to.include <- names(y[i])
        Gene.names.for.table <- which(gene.to.include == MAGMA.gene.regions.for.chromosome$Gene)
        SNPs.to.include <- y[[i]][2:length(y[[i]])]
        
        for (f in 1:length(unrecorded_SNPs)) {
          Base.pairs.index <- which(clumped_SNPs$SNP == SNPs.to.include[unrecorded_SNPs[f]])
          new.row <- list(MAGMA.gene.regions.for.chromosome$CHR[Gene.names.for.table], 
                          SNPs.to.include[unrecorded_SNPs[f]],
                          clumped_SNPs$GD[Base.pairs.index],
                          clumped_SNPs$BP[Base.pairs.index],
                          clumped_SNPs$A1[Base.pairs.index],
                          clumped_SNPs$A2[Base.pairs.index],
                          MAGMA.gene.regions.for.chromosome$GENE_NAME[Gene.names.for.table],
                          gene.to.include,
                          MAGMA.gene.regions.for.chromosome$BP_START[Gene.names.for.table],
                          MAGMA.gene.regions.for.chromosome$BP_END[Gene.names.for.table]
          )
          Gene_clumped_SNPs <- rbindlist(list(Gene_clumped_SNPs,new.row))
        }
      } #name of gene not in list then add row of SNPs with gene identifier to the table()
    }
  }
  assign("test_data_frame",Gene_clumped_SNPs,envir = .GlobalEnv)
}
  if (length_of_genes == "extended") {
    for (i in 1:length(y)) {
      # check whether the gene exists in the updated Gene_clumped_SNPs file
      # If it does not, this means that there is a gene inside a gene or the SNPs were only found inside one of the genes 
      index_of_gene <- which(Gene_clumped_SNPs$Gene == names(y[i]))
      if (length(index_of_gene) == 0){
        gene.to.include <- names(y[i])
        # the first elemnt of y[[i]] is the bp start and end so need to take it into account when finding the SNP in question
        SNPs.to.include <- y[[i]][2:length(y[[i]])]
        
        Gene.names.for.table <- which(gene.to.include == MAGMA.gene.regions.for.chromosome$Gene)
        for (f in 1:length(SNPs.to.include)) {
          # clumped_SNPS is the original unaltered Bim file for each chromosome and pathway
          Base.pairs.index <- which(clumped_SNPs$SNP == SNPs.to.include[f])
          # add a new row corresponding to the Gene_clumped_SNPs format
          new.row <- list(MAGMA.gene.regions.for.chromosome$CHR[Gene.names.for.table], 
                          SNPs.to.include[f], 
                          clumped_SNPs$GD[Base.pairs.index],
                          clumped_SNPs$BP[Base.pairs.index],
                          clumped_SNPs$A1[Base.pairs.index],
                          clumped_SNPs$A2[Base.pairs.index],
                          MAGMA.gene.regions.for.chromosome$GENE_NAME[Gene.names.for.table],
                          gene.to.include,
                          MAGMA.gene.regions.for.chromosome$BP_START[Gene.names.for.table],
                          MAGMA.gene.regions.for.chromosome$BP_END[Gene.names.for.table],
                          MAGMA.gene.regions.for.chromosome$GENE_NAME[Gene.names.for.table],
                          gene.to.include,
                          MAGMA.gene.regions.for.chromosome$BP_start_extended[Gene.names.for.table],
                          MAGMA.gene.regions.for.chromosome$BP_end_extended[Gene.names.for.table])
          Gene_clumped_SNPs <- rbindlist(list(Gene_clumped_SNPs,new.row))
        }
      }else{
        # check if some SNPs are missing
        gene.names.in.y <- Gene_clumped_SNPs[index_of_gene]
        unrecorded_SNPs <- which(!y[[i]][2:length(y[[i]])] %in% gene.names.in.y$SNP)
        # if they are, identify them and do the same for if the gene was missing
        if (length(unrecorded_SNPs) != 0 ) {
          gene.to.include <- names(y[i])
          Gene.names.for.table <- which(gene.to.include == MAGMA.gene.regions.for.chromosome$Gene)
          SNPs.to.include <- y[[i]][2:length(y[[i]])]
          
          for (f in 1:length(unrecorded_SNPs)) {
            Base.pairs.index <- which(clumped_SNPs$SNP == SNPs.to.include[unrecorded_SNPs[f]])
            new.row <- list(MAGMA.gene.regions.for.chromosome$CHR[Gene.names.for.table], 
                            SNPs.to.include[unrecorded_SNPs[f]],
                            clumped_SNPs$GD[Base.pairs.index],
                            clumped_SNPs$BP[Base.pairs.index],
                            clumped_SNPs$A1[Base.pairs.index],
                            clumped_SNPs$A2[Base.pairs.index],
                            MAGMA.gene.regions.for.chromosome$GENE_NAME[Gene.names.for.table],
                            gene.to.include,
                            MAGMA.gene.regions.for.chromosome$BP_START[Gene.names.for.table],
                            MAGMA.gene.regions.for.chromosome$BP_END[Gene.names.for.table],
                            MAGMA.gene.regions.for.chromosome$GENE_NAME[Gene.names.for.table],
                            gene.to.include,
                            MAGMA.gene.regions.for.chromosome$BP_start_extended[Gene.names.for.table],
                            MAGMA.gene.regions.for.chromosome$BP_end_extended[Gene.names.for.table])
            Gene_clumped_SNPs <- rbindlist(list(Gene_clumped_SNPs,new.row))
          }
        } #name of gene not in list then add row of SNPs with gene identifier to the table()
      }
    }
    assign("test_data_frame_extended", Gene_clumped_SNPs,envir = .GlobalEnv)
  }  
}


## Read in maximum Gene limits for each chromosome
# UCSC_file <- "~/Dropbox/whole_genome_testing/output/UCSC_hg19_chromeinfo_length_of_chromosomes.txt"
UCSC_chromosome_lengths <- fread(UCSC_file)
setnames(UCSC_chromosome_lengths,c("CHR","BP","hyperlink"))
UCSC_chromosome_lengths <- UCSC_chromosome_lengths[,.(CHR,BP)]
  
## Sex chromosomes involved?
# if (argument_for_sex_chromosomes == F) {
UCSC_chromosome_lengths <- UCSC_chromosome_lengths[!(CHR == "chrX" | CHR == "chrY")] 
UCSC_chromosome_lengths <- UCSC_chromosome_lengths[mixedorder(UCSC_chromosome_lengths$CHR)]
#}
  
  
## Read in clumped data
# Validation_dataset <- "~/Dropbox/whole_genome_testing/output/CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS.bim"
clumped_SNPs <- fread(Validation_dataset)
setnames(clumped_SNPs,c("CHR","SNP", "GD", "BP", "A1", "A2"))
wd <- getwd()
  

## read in MAGMA's the gene regions 
# NCBI_file <- "~/Dropbox/whole_genome_testing/output/NCBI37.3.gene.loc"
MAGMA.gene.regions <- fread(NCBI_file,colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
setnames(MAGMA.gene.regions, c("Gene","CHR","BP_START","BP_END","STRAND","GENE_NAME"))
  
current_table_name <- MAGMA.gene.regions
# equivalent to MAGMAs window option without ignoring strand
  
current_table_name <- current_table_name[STRAND == "+", BP_start_extended := BP_START - 35000]
current_table_name <- current_table_name[STRAND == "+", BP_end_extended := BP_END + 10000]
current_table_name <- current_table_name[STRAND == "-", BP_start_extended := BP_START - 10000]
current_table_name <- current_table_name[STRAND == "-", BP_end_extended := BP_END + 35000]
  
if (any(current_table_name$BP_start_extended < 0)) {
  reset.to.zero <- which(current_table_name$BP_start_extended < 0)
  current_table_name <- current_table_name[reset.to.zero,BP_start_extended := 0]
}
  
for (ichr in 1:nrow(UCSC_chromosome_lengths)) {
  current_table_chromosome_specific <- current_table_name[CHR == ichr]
    
  if (any(current_table_chromosome_specific$BP_end_extended > UCSC_chromosome_lengths$BP[ichr])){
    set_chromosome_length_limit <- which (current_table_chromosome_specific$BP_end_extended > UCSC_chromosome_lengths$BP[ichr])
    current_table_chromosome_specific <- current_table_chromosome_specific[set_chromosome_length_limit, UCSC_chromosome_lengths$BP[ichr]]
    assign(x = paste0("MAGMA.gene.regions.chromosome", ichr), value = current_table_chromosome_specific, envir = .GlobalEnv)
    print(paste0("needed to adjust gene boundaries at chromosome", ichr, " at Genes: ", current_table_chromosome_specific$GENE_NAME[set_chromosome_length_limit]))
  }else{
    assign(x = paste0("MAGMA.gene.regions.chromosome", ichr), value = current_table_chromosome_specific, envir = .GlobalEnv)
  }
} 

  
  
# Now start on the chromosomes
for (l in chromosomes_to_analyse){
    
  temp_whole_genome_table <- clumped_SNPs[CHR == l]
  temp_gene_annotation_table <- eval(parse(text = paste0("MAGMA.gene.regions.chromosome",l)))
    
  if(nrow(temp_whole_genome_table) == 0) {
    message <- paste0("The number of SNPs on chromosome ", l, "is empty")
    write(message, file = "./output/whole_genome_analysis_info.txt", append = T)
    next()
  }
    
  # Find clumped SNPs
  BP.clumped.SNPs <- temp_whole_genome_table$BP
  
  # Output depending on the arguments
  outputfilename <- paste0("testing_gene_output_chromosome_",l)
  
  Assigning_genes(MAGMA_file = temp_gene_annotation_table, clumped_SNPs = temp_whole_genome_table, BP.clumped.SNPs = BP.clumped.SNPs, outputfilename = outputfilename, gene.regions = "both", chromosome.number = l)
}
  
  
# Read in the 22 tables here and continue!
###
#  
#  for (i in 1:22){
#    assign(paste0("normal_gene_regionstesting_gene_output_chromosome_", i),fread(paste0("./output/normal_gene_regionstesting_gene_output_chromosome_", i, ".txt")), envir = .GlobalEnv)
#  }

# Collect all the tables together into one for the whole genome
# Specify if you wish to include the extended regulatory regions as well

if(Gene_regions_to_analyse == "normal" | Gene_regions_to_analyse == "both"){
  genes_to_snps_list = list()
  for (i in chromosomes_to_analyse) {
    genes_to_snps_list[[i]] <- eval(parse(text = paste0("e$Gene_clumped_SNPs_chromosome", i)))
  }
  Gene_clumped_SNPs <- rbindlist(genes_to_snps_list)
}

if(Gene_regions_to_analyse == "expanded" | Gene_regions_to_analyse == "both"){
 genes_to_snps_list = list()
 for (i in chromosomes_to_analyse) {
   genes_to_snps_list[[i]] <- eval(parse(text = paste0("e$Expanded_Gene_clumped_SNPs_chromosome_", i)))
 }
 Expanded_Gene_clumped_SNPs <- rbindlist(genes_to_snps_list)
 }

  
# read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
    # GENES_to_snps <- scan(file = paste0("./output/CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS.genes.annot"), what = "", sep = "\n")
if(Gene_regions_to_analyse == "normal" | Gene_regions_to_analyse == "both"){
  
    GENES_to_snps <- scan(file = paste0(Validation_dataset,".genes.annot"), what = "", sep = "\n")
    y <- strsplit(GENES_to_snps, "[[:space:]]+")
    names(y) <- sapply(y, '[[', 1)
    y <- lapply(y, '[', -1)
    
    y[[1]] <- NULL
    y[[1]] <- NULL
    
    adding_unread_genes(Gene_clumped_SNPs = Gene_clumped_SNPs, MAGMA.gene.regions.for.chromosome = MAGMA.gene.regions, clumped_SNPs = clumped_SNPs, y = y, chromosome.number = l, length_of_genes = "normal")  
    
    names(test_data_frame) <- c("CHR" ,"SNP", "GD", "BP", "A1", "A2", "Gene_name","GENE", "BP_START","BP_END")
    setcolorder(test_data_frame, c("CHR","SNP","BP","GENE","BP_START","BP_END","GD","A1","A2","Gene_name"))
    
    which(duplicated(test_data_frame$SNP,fromLast = T))
    
    write.table(test_data_frame[, c(1:6), with = F], file = paste0("./output/MAGMA_Gene_regions_for_",Validation_name, "_", Training_name,"txt"), quote = F, row.names = F, col.names = F)
    write(unique(test_data_frame$Gene_name),file = paste0("./output/",Training_name, "_", Validation_name, "_unique_genes.txt"))
    

}

if (Gene_regions_to_analyse == "expanded" | Gene_regions_to_analyse == "both"){
  
  GENES_to_snps <- scan(file = paste0(Validation_dataset,"_extended.genes.annot"), what = "", sep = "\n")
  y_extended <- strsplit(GENES_to_snps, "[[:space:]]+")
  names(y_extended) <- sapply(y_extended, '[[', 1)
  y_extended <- lapply(y_extended, '[', -1)
  
  y_extended[[1]] <- NULL
  y_extended[[1]] <- NULL
  
  adding_unread_genes(Gene_clumped_SNPs = Expanded_Gene_clumped_SNPs, MAGMA.gene.regions.for.chromosome = MAGMA.gene.regions, clumped_SNPs = clumped_SNPs, y = y_extended, chromosome.number = l, length_of_genes = "extended")
  
  test_data_frame_extended <- test_data_frame_extended[,.(CHR, SNP, GD, BP, A1, A2, GENE_NAME, Gene, BP_start_extended, BP_end_extended)]
  names(test_data_frame_extended) <- c("CHR" ,"SNP", "GD", "BP", "A1", "A2", "Gene_name", "GENE", "BP_START", "BP_END")
  setcolorder(test_data_frame_extended, c("CHR","SNP","BP","GENE","BP_START","BP_END","GD","A1","A2","Gene_name"))
  
  which(duplicated(test_data_frame_extended$SNP,fromLast = T))
  
  write.table(test_data_frame_extended[, c(1:6), with = F], file = paste0("./output/MAGMA_Gene_regions_for_",Validation_name, "_", Training_name,"_extended.txt"), quote = F, row.names = F, col.names = F)
  write(unique(test_data_frame_extended$Gene_name),file = paste0("./output/",Training_name, "_", Validation_name, "_unique_genes_extended.txt"))
}

