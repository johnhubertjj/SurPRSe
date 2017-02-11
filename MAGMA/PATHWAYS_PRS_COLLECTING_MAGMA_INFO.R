# Reading_in_Pathway_files
library(data.table)
### environment for functions
e <- new.env()

### Function which assigns genes to the SNP data ####
Assigning_genes <- function(pathway_input, BP.clumped.SNPs, clumped_SNPs, outputfilename, gene.regions = c("normal", "extended")){
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
    browser()
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

##### Setting_up_multi-platform ####
system_information<-Sys.info()

if (system_information[1] == "Windows") fpath <-  "/Users/JJ/" else fpath <-"/Users/johnhubert/"
#####
pathway_sets <- fread(paste0(fpath,"Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/Pocklington2015_134sets_LoFi.txt"))
pathway_sets2 <- fread(paste0(fpath, "Documents/PhD_scripts/Schizophrenia_PRS_pipeline_scripts/GeneWide_BGS_strength.txt"))
useful_pathways <- c("FMRP_targets", "abnormal_behavior", "abnormal_nervous_system_electrophysiology", "abnormal_learning|memory|conditioning", "abnormal_CNS_synaptic_transmission", "Cav2_channels", "abnormal_synaptic_transmission", "5HT_2C", "abnormal_long_term_potentiation", "abnormal_motor_capabilities|coordination|movement", "abnormal_behavioral_response_to_xenobiotic", "abnormal_associative_learning", "Lek2015_LoFintolerant_90", "BGS_top2_mean", "BGS_top2_max")
setnames(pathway_sets, c("V1", "V2"), c("Pathway", "Gene"))
setnames(pathway_sets2, c("V1", "V2"), c("Pathway", "Gene"))
pathway_sets <- merge(pathway_sets, pathway_sets2, by = c("Pathway","Gene"), all = T)
setkey(pathway_sets, Pathway)

for (i in 1:length(useful_pathways)) {
  assign(useful_pathways[i], subset(pathway_sets, Pathway == useful_pathways[i]), envir = .GlobalEnv)
} 

MAGMA.gene.regions <- fread("NCBI37.3.gene.loc",colClasses = c("numeric","character",rep("numeric",2),rep("character",2)))
setnames(MAGMA.gene.regions, c("Gene","CHR","BP_START","BP_END","STRAND","GENE_NAME"))

for (i in 1:length(useful_pathways)) {
  assign(paste0("merged",useful_pathways[i]), merge(eval(parse(text = paste0("`",useful_pathways[i],"`"))), MAGMA.gene.regions, by = "Gene", all = F, sort = F), envir = .GlobalEnv)
  setkey(eval(parse(text = paste0("`","merged",useful_pathways[i],"`"))),STRAND) 
  current_table_name <- eval(parse(text = paste0("`","merged",useful_pathways[i],"`")))
  current_table_name <- current_table_name[STRAND == "+",BP_start_extended := BP_START - 35000]
  current_table_name <- current_table_name[STRAND == "+",BP_end_extended := BP_END + 10000]
  current_table_name <- current_table_name[STRAND == "-",BP_start_extended := BP_START - 10000]
  current_table_name <- current_table_name[STRAND == "-",BP_end_extended := BP_END + 35000]
  setkey(current_table_name, CHR)
  current_table_name <- current_table_name[!"X"]
  assign(paste0("Gene_regions_all_",useful_pathways[i]), current_table_name, envir = .GlobalEnv)
  
for (l in 1:22){
  selecting_chromosomes <- fread(paste0(fpath,"Dropbox/testing_PRS_chromosome_22/output/CLOZUK_GWAS_BGE_chr22_magma_input.bim"))
  names(selecting_chromosomes) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
  temp_pathway_table <- current_table_name[CHR == l]
  selecting_chromosomes_BP <- selecting_chromosomes$BP
  Assigning_genes(pathway_input = temp_pathway_table, clumped_SNPs = selecting_chromosomes, BP.clumped.SNPs = selecting_chromosomes_BP, outputfilename = "testing_gene_output", chromosome.number = l, gene.regions = "extended")
}
}  
selecting_chromosomes <- fread(paste0(fpath,"Dropbox/testing_PRS_chromosome_22/output/CLOZUK_GWAS_BGE_chr22_magma_input.bim"))
