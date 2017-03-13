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

# Set up library
library("data.table")
library("gtools")

# Set working directory
setwd(".")

###############################
## Location Finding Function ##
###############################

Batch.assignment <- function(wherami) {
  if (whereami == "johnhubert") {
    assign("chromosome.number", 22, envir = .GlobalEnv)
    
  } else if (whereami == 'JJ') {
    assign("chromosome.number", 22, envir = .GlobalEnv)
    
  } else if (whereami == "c1020109") {
    
    # Preparing to run in Job array
    AI <- Sys.getenv("PBS_ARRAY_INDEX")
    AI <- as.numeric(AI)
    assign("chromosome.number", AI, envir = .GlobalEnv)
    
  } else {
    stop("current environment NOT at work/home or on servers, please add to script above to specify where you are and what type of analysis you want to do")
  }
}

# Find out whether you want to run the script in a batch format or in serial

### environment for functions
e <- new.env()

###################################################
## Function which assigns genes to the SNP data ##
###################################################
Assigning_genes <- function(MAGMA_file, BP.clumped.SNPs, clumped_SNPs, outputfilename, gene.regions = c("normal", "extended"), chromosome.number = l){
  
  if(gene.regions == "normal"){
    
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
                file = paste0("./output/normal_gene_regions", outputfilename, ".txt"), 
                quote = F, 
                row.names = F)
    assign(paste0("Gene_clumped_SNPs_chromosome",l), Gene_clumped_SNPs, envir = e)
  }
  
  if(gene.regions == "extended"){
    
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
                file = paste0("./output/", outputfilename, ".txt"), 
                quote = F, 
                row.names = F)
    
    assign("Expanded_Gene_clumped_SNPs", Gene_clumped_SNPs, envir = e)
  }
}

###############################################################################
### Function which adds duplicate SNPs which happen to be inside other genes ##
###############################################################################
adding_unread_genes <- function(Gene_clumped_SNPs, MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number, args = c('Batch', 'Serial')){

  for (i in 1:length(y)) {
    index_of_gene <- which(Gene_clumped_SNPs$Gene_ID == names(y[i]))
    if (length(index_of_gene) == 0){
      gene.to.include <- names(y[i])
      SNPs.to.include <- y[[i]][2:length(y[[i]])]
      
      Gene.names.for.table <- which(gene.to.include == MAGMA.gene.regions.for.chromosome$ID)
      for (f in 1:length(SNPs.to.include)) {
        Base.pairs.index <- which(clumped_SNPs$SNP == SNPs.to.include[f])

        if (args == 'Serial') {
          chromosome.number <- sub("(?<=^d+)\\:.*","",SNPs.to.include, perl = F) # Not sure if this works, test it 
        } 
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

# Run script if you are running as a Job array #

if( args[1] == 'Batch') {
  # Assign chromosome number as job array variable
  Batch.assignment(whereami)


  ## Read in clumped data
  clumped_SNPs <- fread("combined_CLOZUK_PGC_CLUMPED_FINAL.txt")
  wd <-getwd()

  ## read in MAGMA's the gene regions 
  MAGMA.gene.regions <- fread("NCBI37.3.gene.loc",colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
  colnames(MAGMA.gene.regions) <- c("ID","Chromosome","BP_1","BP_2","strand","Gene_symbol")
  
  ## Limit the data down to the specific chromsome the array is on
  Index.for.chromosome <- MAGMA.gene.regions[,.I[grep(paste0("^",chromosome.number),Chromosome, perl = T, invert = F)]]
  MAGMA.gene.regions.for.chromosome <- MAGMA.gene.regions[Index.for.chromosome]
  MAGMA.gene.regions.for.chromosome <- MAGMA.gene.regions.for.chromosome[,c('ID','Chromosome','BP_1','BP_2','Gene_symbol'), with = F]
  
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

  write.table(test_data_frame[, c(1:6), with = F], file = paste0("./output/MAGMA_Gene_regions_for_python_script_chr_",chromosome.number,".txt"), quote = F, row.names = F)
  write(unique(test_data_frame$Gene_name), file = paste0("./output/PGC_CLOZUK_unique_genes_chr_",chromosome.number,".txt"))
}

# Run script if you are running in serial #
if (args[1] == 'Serial' ) {
  
  ## Read in maximum Gene limits for each chromosome
  UCSC_chromosome_lengths <- fread("./output/UCSC_hg19_chromeinfo_length_of_chromosomes.txt")
  setnames(UCSC_chromosome_lengths,c("CHR","BP","hyperlink"))
  UCSC_chromosome_lengths <- UCSC_chromosome_lengths[,.(CHR,BP)]
  
  ## Sex chromosomes involved?
  # if (argument_for_sex_chromosomes == F) {
  UCSC_chromosome_lengths <- UCSC_chromosome_lengths[!(CHR == "chrX" | CHR == "chrY")] 
  UCSC_chromosome_lengths <- UCSC_chromosome_lengths[mixedorder(UCSC_chromosome_lengths$CHR)]
  #}
  
  
  ## Read in clumped data
  clumped_SNPs <- fread("./output/CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPs.bim")
  setnames(clumped_SNPs,c("CHR","SNP", "GD", "BP", "A1", "A2"))
  wd <- getwd()
  

  ## read in MAGMA's the gene regions 
  MAGMA.gene.regions <- fread("./output/NCBI37.3.gene.loc",colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
  setnames(MAGMA.gene.regions, c("Gene","CHR","BP_START","BP_END","STRAND","GENE_NAME"))
  
  current_table_name <- MAGMA.gene.regions
  # equivalent to MAGMAs window option without ignoring strand
  
  current_table_name <- current_table_name[STRAND == "+", BP_start_extended := BP_START - 35000]
  current_table_name <- current_table_name[STRAND == "+", BP_end_extended := BP_END + 10000]
  current_table_name <- current_table_name[STRAND == "-", BP_start_extended := BP_START - 10000]
  current_table_name <- current_table_name[STRAND == "-",BP_end_extended := BP_END + 35000]
  
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
  for (l in 1:22){
    
    temp_whole_genome_table <- clumped_SNPs[CHR == l]
    temp_gene_annotation_table <- eval(parse(text = paste0("MAGMA.gene.regions.chromosome",l)))
    
    if(nrow(temp_whole_genome_table) == 0) {
      message <- paste0("The number of SNPs on chromosome ",l,"is empty")
      write(message, file = "whole_genome_analysis_info.txt", append = T)
      next()
    }
    
    # Find clumped SNPs
    BP.clumped.SNPs <- temp_whole_genome_table$BP
    outputfilename_normal <- paste0("testing_gene_output_chromosome_",l)
    Assigning_genes(MAGMA_file = temp_gene_annotation_table, clumped_SNPs = temp_whole_genome_table, BP.clumped.SNPs = BP.clumped.SNPs, outputfilename = outputfilename_normal, gene.regions = "normal")
  }
  
  
# Read in the 22 tables here and continue!
###
  
    # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
    GENES_to_snps <- scan(file = paste0("./output/CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS.genes.annot"), what = "", sep = "\n")
    y <- strsplit(GENES_to_snps, "[[:space:]]+")
    names(y) <- sapply(y, '[[', 1)
    y <- lapply(y, '[', -1)
    
    y[[1]] <- NULL
    y[[1]] <- NULL
    
    adding_unread_genes(Gene_clumped_SNPs = e$Expanded_Gene_clumped_SNPs, MAGMA.gene.regions.for.chromosome = temp_pathway_table, clumped_SNPs = selecting_chromosomes, y = y, chromosome.number = l, pathway = useful_pathways[i])  
    SNPs_for_clumping <- paste0("chromosome_",l,"_", useful_pathways[i],"SNPs_for_clumping")
    assign(SNPs_for_clumping, unique(test_data_frame$SNP), envir = .GlobalEnv)
    write.table(test_data_frame,file = paste0("./output/",useful_pathways[i],"_chromosome_",l,"_extended_data_table.txt"),quote = F, row.names = F)
    SNPs_for_clumping <- paste0("`","chromosome_",l,"_",useful_pathways[i],"SNPs_for_clumping","`")
    write(eval(parse(text = SNPs_for_clumping)), file = paste0("./output/chromosome_",l,"_", useful_pathways[i],"SNPs_for_clumping.txt"))
    rm(list = paste0("chromosome_",l,"_", useful_pathways[i],"SNPs_for_clumping"))
  }
}

  
  # Create table and write to file for reference
  Assigning_genes(MAGMA.gene.regions, 
                  BP.clumped.SNPs = BP.clumped.SNPs, 
                  clumped_SNPs = clumped_SNPs, 
                  outputfilename = "Genomic_1000kb_r2_zeropoint2_PGC_CLOZUK")
  
  # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
  GENES_to_snps <- scan(file = paste0("./output/CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS.genes.annot"), what = "", sep = "\n")
  y <- strsplit(GENES_to_snps, "[[:space:]]+")
  names(y) <- sapply(y, '[[', 1)
  y <- lapply(y, '[', -1)
  
  y[[1]] <- NULL
  y[[1]] <- NULL
  
  chromosome.number <- "NA"
  
  adding_unread_genes(e$test_clumped_snps,MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number)  
  
  names(test_data_frame) <- c("CHR" ,"SNP", "BP", "P", "Gene_name","BP_START","BP_END","GENE")
  setcolorder(test_data_frame, c("CHR","SNP","BP","GENE","BP_START","BP_END","P","Gene_name"))
  
  which(duplicated(test_data_frame$SNP,fromLast = T))
  
  write.table(test_data_frame[, c(1:6), with = F], file = paste0("./output/MAGMA_Gene_regions_for_python_script_chr_",chromosome.number,".txt"), quote = F, row.names = F)
  write(unique(test_data_frame$Gene_name), file = paste0("./output/PGC_CLOZUK_unique_genes_chr_",chromosome.number,".txt"))
}
########################################
########## Pathway Results #############
########################################

# Reading_in_Pathway_files
library(data.table)
### environment for functions
e <- new.env()

### Function which assigns genes to the SNP data ####
### Currently designed for the input from bim files ###
Assigning_genes <- function(pathway_input, BP.clumped.SNPs, clumped_SNPs, outputfilename, gene.regions = c("normal", "extended"),chromosome.number = l){
  GR <- deparse(substitute(pathway_input))
  
  if(gene.regions == "MAGMA.gene.regions.for.chromosome"){
    
    test_vector <- rep(0,nrow(clumped_SNPs))
    gene_ID_vector <- rep(0,nrow(clumped_SNPs))
    BP1 <- rep(0,nrow(clumped_SNPs))
    BP2 <- rep(0,nrow(clumped_SNPs))
    
    for (i in 1:nrow(clumped_SNPs)) {
      tmp_gene <- which ((pathway_input$BP_START <= BP.clumped.SNPs[i] & pathway_input$BP_START >= BP.clumped.SNPs[i]) == T)
      if (length(tmp_gene) == 1) {
        test_vector[i] <- pathway_input$GENE_NAME[tmp_gene]
        gene_ID_vector[i] <- pathway_input$Gene[tmp_gene]
        BP1 <- pathway_input$BP1[tmp_gene]
        BP2 <- pathway_input$BP2[tmp_gene]
        
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
  
  if(gene.regions == "extended"){
    
    test_vector <- rep(0,nrow(clumped_SNPs))
    gene_ID <- rep(0,nrow(clumped_SNPs))
    BP1 <- rep(0,nrow(clumped_SNPs))
    BP2 <- rep(0,nrow(clumped_SNPs))
    
    for (i in 1:nrow(clumped_SNPs)) {
      tmp_gene <- which ((pathway_input$BP_start_extended <= BP.clumped.SNPs[i] & pathway_input$BP_end_extended >= BP.clumped.SNPs[i]) == T)
      
      if (length(tmp_gene) == 1) {
        test_vector[i] <- pathway_input$GENE_NAME[tmp_gene]
        gene_ID[i] <- pathway_input$Gene[tmp_gene]
        BP1[i] <- pathway_input$BP_start_extended[tmp_gene]
        BP2[i] <- pathway_input$BP_end_extended[tmp_gene]
        
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
                file = paste0("./output/", outputfilename, chromosome.number, ".txt"), 
                quote = F, 
                row.names = F)
    
    assign("Expanded_Gene_clumped_SNPs", Gene_clumped_SNPs, envir = e)
  }
}

### Function which adds duplicate SNPs which happen to be inside other genes
### Currently designed for extended gene regions ONLY
### requires the run of gene.out R script to obtain the specific gene.annot file for each pathway and chromosome
adding_unread_genes <- function(Gene_clumped_SNPs, MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number, pathway){
  # y is the input gene list from MAGMA
  # Gene_clumped_SNPs is the updated bim table from the initial assignment of SNPs
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
        new.row <- list(chromosome.number, 
                        SNPs.to.include[f], 
                        clumped_SNPs$GD[Base.pairs.index],
                        clumped_SNPs$BP[Base.pairs.index],
                        clumped_SNPs$A1[Base.pairs.index],
                        clumped_SNPs$A2[Base.pairs.index],
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
          new.row <- list(chromosome.number, 
                          SNPs.to.include[unrecorded_SNPs[f]],
                          clumped_SNPs$GD[Base.pairs.index],
                          clumped_SNPs$BP[Base.pairs.index],
                          clumped_SNPs$A1[Base.pairs.index],
                          clumped_SNPs$A2[Base.pairs.index],
                          MAGMA.gene.regions.for.chromosome$GENE_NAME[Gene.names.for.table],
                          gene.to.include,
                          MAGMA.gene.regions.for.chromosome$BP_start_extended[Gene.names.for.table],
                          MAGMA.gene.regions.for.chromosome$BP_end_extended[Gene.names.for.table]
          )
          Gene_clumped_SNPs <- rbindlist(list(Gene_clumped_SNPs,new.row))
        }
      } #name of gene not in list then add row of SNPs with gene identifier to the table()
    }
  }
  assign("test_data_frame",Gene_clumped_SNPs,envir = .GlobalEnv)
}

##### Setting_up_multi-platform ####
system_information<-Sys.info()

if (system_information[1] == "Windows") fpath <-  "/Users/JJ/" else fpath <-"/Users/johnhubert/"

##### read in pathway sets from Antonio
pathway_sets <- fread(paste0(fpath,"Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Pocklington2015_134sets_LoFi.txt"))
pathway_sets2 <- fread(paste0(fpath, "Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/GeneWide_BGS_strength.txt"))

# input the name of the pathways from the orignal file
useful_pathways <- c("FMRP_targets", "abnormal_behavior", "abnormal_nervous_system_electrophysiology", "abnormal_learning|memory|conditioning", "abnormal_CNS_synaptic_transmission", "Cav2_channels", "abnormal_synaptic_transmission", "5HT_2C", "abnormal_long_term_potentiation", "abnormal_motor_capabilities|coordination|movement", "abnormal_behavioral_response_to_xenobiotic", "abnormal_associative_learning", "Lek2015_LoFintolerant_90", "BGS_top2_mean", "BGS_top2_max")
setnames(pathway_sets, c("V1", "V2"), c("Pathway", "Gene"))
setnames(pathway_sets2, c("V1", "V2"), c("Pathway", "Gene"))

# merge the two input pathway files together
pathway_sets <- merge(pathway_sets, pathway_sets2, by = c("Pathway","Gene"), all = T)
# set the key to the pathway name
setkey(pathway_sets, Pathway)

# assign an individual table for each pathway
for (i in length(useful_pathways)) {
  assign(useful_pathways[i], subset(pathway_sets, Pathway == useful_pathways[i]), envir = .GlobalEnv)
} 

#read in general MAGMA annotation file
MAGMA.gene.regions <- fread("./output/NCBI37.3.gene.loc",colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
setnames(MAGMA.gene.regions, c("Gene","CHR","BP_START","BP_END","STRAND","GENE_NAME"))


# loop through both chromosomes and pathway tables
for (i in 1:length(useful_pathways)) {
  
  # merge the original MAGMA annotation file and the pathway file in order to get extra info including BP extensions for regulatory regions 
  assign(paste0("merged",useful_pathways[i]), merge(eval(parse(text = paste0("`",useful_pathways[i],"`"))), MAGMA.gene.regions, by = "Gene", all = F, sort = F), envir = .GlobalEnv)
  setkey(eval(parse(text = paste0("`","merged",useful_pathways[i],"`"))),STRAND) 
  current_table_name <- eval(parse(text = paste0("`","merged",useful_pathways[i],"`")))
  
  # equivalent to MAGMAs window option without ignoring strand
  current_table_name <- current_table_name[STRAND == "+", BP_start_extended := BP_START - 35000]
  current_table_name <- current_table_name[STRAND == "+", BP_end_extended := BP_END + 10000]
  current_table_name <- current_table_name[STRAND == "-", BP_start_extended := BP_START - 10000]
  current_table_name <- current_table_name[STRAND == "-", BP_end_extended := BP_END + 35000]
  
  # setkey to the chromosome and remove all SNPs in the X chromosome to save time
  setkey(current_table_name, CHR)
  current_table_name <- current_table_name[!"X"]
  assign(paste0("Gene_regions_all_",useful_pathways[i]), current_table_name, envir = .GlobalEnv)
  
  # Now start on the chromosomes
  for (l in 1:22){
    
    selecting_chromosomes <- fread(paste0(fpath,"Documents/testing_PRS_chromosome_22/test_chr5/output/CLOZUK_GWAS_BGE_chr",l,"_magma_input.bim"))
    names(selecting_chromosomes) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
    temp_pathway_table <- current_table_name[CHR == l]
    if(nrow(temp_pathway_table) == 0) {
      message <- paste0("The number of SNPs for",useful_pathways[i],"on chromosome ",l,"is empty")
      write(message, file = "Pathways_analysis_info.txt", append = T)
      next()
    }
    selecting_chromosomes_BP <- selecting_chromosomes$BP
    Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, outputfilename = "testing_gene_output", chromosome.number = l, gene.regions = "extended")
    
    
    
    # read in MAGMA's input and add any genes which happen to be inside other genes or crossed over with other genes
    GENES_to_snps <- scan(file = paste0("./output/",l,"_CLOZUK_PGC_SNPs_",useful_pathways[i],"pathway.genes.annot"), what = "", sep = "\n")
    y <- strsplit(GENES_to_snps, "[[:space:]]+")
    names(y) <- sapply(y, '[[', 1)
    y <- lapply(y, '[', -1)
    
    y[[1]] <- NULL
    y[[1]] <- NULL
    
    adding_unread_genes(Gene_clumped_SNPs = e$Expanded_Gene_clumped_SNPs, MAGMA.gene.regions.for.chromosome = temp_pathway_table, clumped_SNPs = selecting_chromosomes, y = y, chromosome.number = l, pathway = useful_pathways[i])  
    SNPs_for_clumping <- paste0("chromosome_",l,"_", useful_pathways[i],"SNPs_for_clumping")
    assign(SNPs_for_clumping, unique(test_data_frame$SNP), envir = .GlobalEnv)
    write.table(test_data_frame,file = paste0("./output/",useful_pathways[i],"_chromosome_",l,"_extended_data_table.txt"),quote = F, row.names = F)
    SNPs_for_clumping <- paste0("`","chromosome_",l,"_",useful_pathways[i],"SNPs_for_clumping","`")
    write(eval(parse(text = SNPs_for_clumping)), file = paste0("./output/chromosome_",l,"_", useful_pathways[i],"SNPs_for_clumping.txt"))
    rm(list = paste0("chromosome_",l,"_", useful_pathways[i],"SNPs_for_clumping"))
  }
}



names(test_data_frame) <- c("CHR" ,"SNP", "BP", "P", "Gene_name","BP_START","BP_END","GENE")
setcolorder(test_data_frame, c("CHR","SNP","BP","GENE","BP_START","BP_END","P","Gene_name"))

which(duplicated(test_data_frame$SNP,fromLast = T))

write.table(test_data_frame[, c(1:6), with = F], file = paste0("./output/MAGMA_Gene_regions_for_python_script_chr_",chromosome.number,".txt"), quote = F, row.names = F)
write(unique(test_data_frame$Gene_name),file = paste0("./output/PGC_CLOZUK_unique_genes_chr_",chromosome.number,".txt"))




