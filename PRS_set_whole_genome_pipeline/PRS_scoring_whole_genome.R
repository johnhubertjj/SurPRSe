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

# Specify the different input tables #
Training_name <- args[3]
Validation_name <- args [4]
bim_file_path <- args[5]
ending_name <- args[6]

# This will be a problem, write in a script that defines what the arguments are (probably have to be linked with the arguments script however...stop unneeded analysis)
significance_thresholds <- as.numeric(args[c(7:length(args))])
print(significance_thresholds)

# Read in data
PGC <- fread(paste0("./", Training_name, "_", Validation_name,"_output/combined_", Training_name, "_table_with_CHR.POS_identifiers.txt"))
test_data_frame <- fread(bim_file_path)
setnames(test_data_frame, c("CHR", "SNP","Gene_ID","BP","A1", "A2"))

# rename (because it was already here, plus I want to see the both the reaction of programmers or whether they notice when I re-assign memo\iry for no reason)
p.value.thresholds <- significance_thresholds

## check that it is merging properly here (after analysis is run)
combined.CLOZUK.PGC.clumped.Genomic.SNPs <- merge(test_data_frame, PGC, by.x="SNP", by.y="SNP", all=F, sort=F)
combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1 <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1.y)
combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2 <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2.y)

# Produce a score file for each p-value threshold and write to the new PRS_scoring file
for (i in 1:length(p.value.thresholds)) {
  
    a <- copy(combined.CLOZUK.PGC.clumped.Genomic.SNPs)    
    SNPs <- a[, .I[which(P <= p.value.thresholds[i])]]
    
    if (length(SNPs) != 0){
      a <- a[SNPs, .(SNP,A1.y,BETA)]
      
      filename <- paste0("./", Training_name, "_", Validation_name,"_output/PRS_scoring/", Training_name, "_", Validation_name, "_", ending_name, "_significance_threshold_at_", p.value.thresholds[i],".score")
      
      write.table(file = filename, a, row.names = F, col.names = F, quote = F, sep="\t")
}
}  
    #command <- paste0("plink --bfile ./", Training_name, "_", Validation_name,"_output/",Validation_name,"_",Training_name, "_FULL_GENOME_CLUMPED --score ", filename, " --out ./", Training_name, "_", Validation_name,"_output/PRS_scoring/",Training_name,"_",Validation_name,"_whole_genome_significance_threshold_at_",p.value.thresholds[i])
      #system(command)
      #rm(a)
    #}else{
      #next()
    #}
#}

