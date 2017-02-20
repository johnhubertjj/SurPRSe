### gain the number of files that match the chromosomes used ###
setwd("/Users/johnhubert/Dropbox/PGC/CLOZUK_PGC_raven_output/")
current_files <- system("ls",intern = T)
current_files <- c(current_files,1,"ahroghouhguo[1]")
current_files_index <- grep("424000\\[\\d+\\].raven0.OU",current_files,perl = T)
current_files <- current_files[current_files_index]

### For loop (change to apply to batch script "PBS_NUM_ARRAY") ####
file_collection <- list()
for (i in 1:length(current_files)) {
  file_collection[[i]] <- read.delim(paste0("424000[",i,"].raven0.OU"))
} 