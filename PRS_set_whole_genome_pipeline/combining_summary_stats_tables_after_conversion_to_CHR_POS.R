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
Chromosomes_to_split <- as.numeric(args[c(5:length(args))])
print(Chromosomes_to_split)

## Read in the latest Summary stats tables after converting to one table
for (i in Chromosomes_to_split){
  assign(paste0(Training_name, "_table", i), fread(paste0("./", Training_name, "_", Validation_name,"_output/", Training_name,"_table",i,"_new.txt")),envir = .GlobalEnv)
}

l = list()

## print out to one table under a common filename
for (i in Chromosomes_to_split) {
  l[[i]] <- eval(parse(text = paste0(Training_name,"_table",i)))
}
combined_final_table <- rbindlist(l)
write.table (combined_final_table, file = paste0("./", Training_name, "_", Validation_name,"_output/combined_", Training_name, "_table_with_CHR.POS_identifiers.txt"), quote = F, row.names = F)

