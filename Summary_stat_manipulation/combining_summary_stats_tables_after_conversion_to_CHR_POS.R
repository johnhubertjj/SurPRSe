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

# specify the different input tables #
Training_name <- args[3]
Is_Training_set_split_up <- args[4]
Chromosomes_to_analyse <- args[5]

## Are the tables split up in some way?
if (Is_Training_set_split_up == "TRUE"){

  ## Create a larger table ##
  for (i in Chromosomes_to_analyse){
  assign(paste0(Training_name,"_table",i),fread(paste0("./output/",Training_name,"_table",i,"_new.txt")),envir = .GlobalEnv)
  }

  l = list()
  
  for (i in 1:22) {
  l[[i]] <- eval(parse(text = paste0(Training_name,"_table",i)))
  }
  combined_final_table <- rbindlist(l)
  
# Write new table to the output directory
  write.table (combined_final_table, file = paste0("./output/combined_",Training_name,"_table_with_CHR.POS_identifiers.txt", quote = F, row.names = F))
}

# End timer
proc.time() - ptm  