# Splitting GeM data
setwd("~/Dropbox/GeM_genotypes/")

GeM <- fread("~/Dropbox/GeM_genotypes/gem.gwa.used.sample.2131.bim")
number_of_chromosomes <- c(1:22)

for (i in number_of_chromosomes){
  current_matrix <- GeM[V1 == i]
  current_SNPs <- current_matrix$V2 
  write(current_SNPs, file = paste0("convert_to_chromosomes_GeM", i, ".txt"))
}