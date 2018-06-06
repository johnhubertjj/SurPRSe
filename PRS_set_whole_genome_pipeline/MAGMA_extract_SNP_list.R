### Start Timer
ptm <- proc.time()

### Library
library(data.table)

## set wd
setwd(".")

########################################
# adding in arguments from BASH script #
########################################
args <- commandArgs(trailingOnly = T)
getwd()
print(args)

## specify the different input tables #
Training_name <- args[3]
Validation_name <- args [4]
chromosome.number <- as.numeric(args[c(5:length(args))])
print(chromosome.number)

## Read in the latest Summary stats tables after converting to one table


for (i in chromosome.number){
  assign(paste0("Training_table", i), fread(paste0("./",Training_name,"_",Validation_name,"_output/",Validation_name,"_chr",i,"_consensus_with_", Training_name,"_flipped_alleles_no_duplicates.bim")))
}

l = list()

## print out to one table under a common filename
for (i in chromosome.number) {
  l[[i]] <- eval(parse(text = paste0("Training_table",i)))
}
combined_final_table <- rbindlist(l)

training_tables_to_merge <- fread(paste0("./", Training_name, "_", Validation_name,"_output/combined_", Training_name, "_table_with_CHR.POS_identifiers.txt"))


colnames(combined_final_table) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
SNP_list_table <- merge(combined_final_table,training_tables_to_merge, by = c("SNP","CHR","BP","A1","A2"), all = F)

cat("The Number of SNPS consistent between both datasets after QC =",nrow(SNP_list_table))

write.table(SNP_list_table, file = paste0( "combined_", Training_name,"_", Validation_name, "_table_with_CHR.POS_identifiers_flipped_alleles_no_duplicates_MAF_INFO.txt"), quote = F, row.names = F, col.names = T)

