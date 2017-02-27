adding_unread_genes <- function(Gene_clumped_SNPs, MAGMA.gene.regions.for.chromosome, clumped_SNPs, y, chromosome.number){
  browser()
  for (i in 1:length(y)) {
    index_of_gene <- which(Gene_clumped_SNPs$Gene_ID == names(y[i]))
    if (i == 10) browser()
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
                        MAGMA.gene.regions.for.chromosome$Gene_symbol[Gene.names.for.table],
                        gene.to.include)
        if (i == 1 & f == 1) {
          rbindlist(list(Gene_clumped_SNPs[1, ], new.row, Gene_clumped_SNPs[i + 1:length(y), ]))
        }else if(i == 1 & f > 1) {
          rbindlist(list(Gene_clumped_SNPs[1 + (f - 1), ], new.row, Gene_clumped_SNPs[(i + f):length(y), ]))
        }else if (i > 1 & f == 1) {
          rbindlist(list(Gene_clumped_SNPs[1:i, ], new.row, Gene_clumped_SNPs[i + 1:length(y), ]))
        }else{
          rbindlist(list(Gene_clumped_SNPs[1:(i + (f - 1)), ], new.row, Gene_clumped_SNPs[(i + f) :length(y), ]))
        }
      }
    }else{
      next()
    } #name of gene not in list then add row of SNPs with gene identifier to the table()
  }
}
