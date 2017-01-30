System_information <- Sys.info()
whereami <- System_information['user']

if (whereami == "johnhubert") {
  Filename <- "/Users/johnhubert/Documents/PGC/daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt"

} else if (whereami == 'JJ') {
  Filename <- "get_file_here!!!!!!!"
  
} else if (whereami == "c1020109") {
  # Preparing to read in arguments from script
  Filename <- Sys.getenv("$1")
  
} else {
  stop("current environment NOT at work/home or on servers, please add to script above to specify where you are and what type of analysis you want to do")
}

## Reading in PGC/CLOZUK script

## prerequisites
library(data.table)

## setting working directory
setwd(".")


## Reading in PGC data
## Select for CHR 22
PGC_data <- fread(Filename)
setkey(PGC_data,CHR)

for (chromosome.number in 1:22) { 
write.table(eval(parse(text = paste0("PGC_data[J(",chromosome.number,")]"))), file = paste0("PGC_table",chromosome.number,".txt"), quote = F, row.names = F)
}

