New_adding_genes_from_magma_function <- function(MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number, pathway){

ptm <- proc.time()  
Table_of_integers <- as.data.frame(matrix(data = NA,nrow = length(y), ncol = 2))
rows_of_tables_to_use <- matrix(data = NA, nrow = length(y), ncol = 2)
rows_of_tables_to_use[1,1] <- 1
#rows_of_tables_to_use[1,2] <- length(y[[i]][2:length(y[[i]])])
current_end <- 0

  for (i in 1:length(y)){
    gene.to.include <- names(y[i])
    SNPs.to.include <- y[[i]][2:length(y[[i]])]
    Table_of_integers[i,1] <- gene.to.include
    Table_of_integers[i,2] <- length(SNPs.to.include)
    length_of_integer_increase <- length(SNPs.to.include)
    if(i == 1){
    current_end <- c(current_end + length_of_integer_increase)
    rows_of_tables_to_use[i,2] <- current_end
    }else{
      current_end <- c(current_end + length_of_integer_increase - 1)
      rows_of_tables_to_use[i,2] <- current_end
    }
    if( i != length(y)){
      current_end <- c(current_end + 1)
      rows_of_tables_to_use[i+1,1] <- current_end
    }
  }

Complete_table <- as.data.table(matrix(data = as.character(NA), nrow = sum(Table_of_integers[,2]), ncol = 10))
 

for (i in 1:length(y)){
  start_integer <- rows_of_tables_to_use[i,1]
  end_integer <- rows_of_tables_to_use[i,2]
  gene.to.include <- names(y[i])
  SNPs.to.include <- y[[i]][2:length(y[[i]])]
  altered_scale_2 <- c(start_integer:end_integer)
  Gene.names.for.table <- which(gene.to.include == MAGMA.gene.regions.for.chromosome$Gene)
  
  for (l in start_integer:end_integer){
    new_index_to_figure_out_selecting_chromosomes_index <- which(altered_scale_2 == l)
    Base.pairs.index <- which(clumped_SNPs$SNP == SNPs.to.include[new_index_to_figure_out_selecting_chromosomes_index])
    
    new.row <- c(chromosome.number, 
                    SNPs.to.include[new_index_to_figure_out_selecting_chromosomes_index], 
                    clumped_SNPs$GD[Base.pairs.index],
                    clumped_SNPs$BP[Base.pairs.index],
                    clumped_SNPs$A1[Base.pairs.index],
                    clumped_SNPs$A2[Base.pairs.index],
                    MAGMA.gene.regions.for.chromosome$GENE_NAME[Gene.names.for.table],
                    gene.to.include,
                    MAGMA.gene.regions.for.chromosome$BP_start_extended[Gene.names.for.table],
                    MAGMA.gene.regions.for.chromosome$BP_end_extended[Gene.names.for.table])
    for (j in seq_len(ncol(Complete_table))){
      set(Complete_table,l,j,new.row[j])
  }
  }
}

#End Timer
proc.time() - ptm

X = list(a=1:5, a=6:10)

adding_unread_genes <- function(MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number){

Melted_MAGMA_list <- melt(y)
names(Melted_MAGMA_list) <- c("SNP", "Gene")
Melted_MAGMA_list$SNP <- as.character(Melted_MAGMA_list$SNP) 
clumped_SNPs$SNP <- as.character(clumped_SNPs$SNP)
merged_table_one <- merge(clumped_SNPs, Melted_MAGMA_list, by = "SNP", all = F)
merged_table_one$Gene <- as.numeric(merged_table_one$Gene)
Gene_clumped_SNPs <- merge(MAGMA.gene.regions,merged_table_one, by = "Gene", all.y = T)
Gene_clumped_SNPs[,"CHR.y" :=NULL]  #remove extra_column
setcolorder(Gene_clumped_SNPs, c("CHR.x","SNP", "GD", "BP", "A1", "A2","GENE_NAME","Gene", "BP_START", "BP_END", "BP_start_extended","BP_end_extended","STRAND"))
setnames(Gene_clumped_SNPs, c("CHR.x","Gene"), c("CHR","GENE_NUMBER"))
assign("test_data_frame", Gene_clumped_SNPs, envir = e)

#End Timer
proc.time() - ptm
}


