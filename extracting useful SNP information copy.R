### Start Timer
ptm <- proc.time()

## Load packages ## 
library(stringr)
library(data.table)
library(plyr)

## Check where you are ## 
System_information <- Sys.info()
whereami <- System_information['user']

if (whereami == "johnhubert") {
  chromosome.number <- args[7]
  
} else if (whereami == 'JJ') {
  chromosome.number <- args[7]
  
} else if (whereami == "c1020109") {
  
  # Preparing to run in Job array
  JOB_ID <- Sys.getenv("JOB_ID")
  chromosome.number <- as.numeric(JOB_ID)
  
} else {
  stop("current environment NOT at work/home or on servers, please add to script above to specify where you are and what type of analysis you want to do")
}

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)

# specify the different input tables #
Training_datatable <- paste0("./output/",args[3],"_new.txt")
Validation_datatable_bim_file <- paste0(args[4],".bim")
Training_name <- args[5]
Validation_name <- args[6]
chromosomes_to_parse <- args[7]
Rout_file_storage <- args[8]
INFO_decision <- args[9]
MAF_decision <- args[10]
MAF_threshold <- args[11]
INFO_decision <- args[12]


### gain the number of files that match the chromosomes used ###
setwd("./" ,Rout_file_storage)
current_files <- system("ls", intern = T)


### For loop (change to apply to batch script "PBS_NUM_ARRAY") ####
file_collection <- list()
for (i in 1:length(chromosomes_to_parse)) {
  current_files_index <- grep(paste0("\\d+\\[",i,"\\].raven0.OU"),current_files,perl = T)
  file_collection[[i]] <- read.delim(current_files[current_files_index])
} 

### Find the rows where cat function appears
informative.data.frame <- list()

for (i in 1:length(chromosomes_to_parse)) {

length_of_combined_datasets_integer <- grep("Chr\\:\\s\\d+",file_collection[[i]][,1],perl = T)
new_vector_two <- as.character(file_collection[[i]][length_of_combined_datasets_integer,1])

### Check for warning messages that have been flagged
warnings <- grep("Warning\\smessage", file_collection[[i]][,1], perl = T)

if(length(warnings) > 0){
warnings <- c(warnings, warnings + 1)
Warnings_here <- as.character(file_collection[[i]][warnings,1])
}

# Extract all strings with a space before a number at the end : should only be a chromosome number for now
abd <-str_extract_all(new_vector_two,"(?<=\\s)\\d+$")

# Print out the number of SNPs throughout the pipeline
if (MAF_decision == "TRUE" & INFO_decision == "TRUE"){
  informative.data.frame[[i]] <- data.frame(info = c(paste0("Number of SNPs in", Training_name, "before isolating MAF > ", MAF_threshold, " and INFO > ", INFO_threshold, " in Chr:", i, " N="),
                                                     paste0("Number of SNPs in", Training_name, "after isolating MAF > ", MAF_threshold, " and INFO > ", INFO_threshold, " in Chr:", i, " N="),
                                                     paste0("Number of SNPS in ", Training_name, "before flipping Chr:", i," N="),
                                                     paste0("Number of SNPs in ", Validation_name, "before flipping Chr:", i," N="),
                                                     paste0("Combined Number of SNPs BEFORE flipping between ",Validation_name," and ",Training_name," Chr:", i," N="),
                                                     paste0(Validation_name,"_", Training_name, " Chr:", i," remove A-T, C-G Step one: N="),
                                                     paste0(Validation_name,"_", Training_name, " Chr:", i," remove A-T, C-G Step two: N="),
                                                     paste0("Combined all A-T,C-G remove and other mismatches Chr", i," N=")),
                                            N = c(abd[[1]], abd[[2]], abd[[3]], abd[[4]], abd[[5]], abd[[6]], abd[[7]], abd[[8]]))
}
if (MAF_decision == "TRUE"){
  informative.data.frame[[i]] <- data.frame(info = c(paste0("Number of SNPs in", Training_name, "before isolating MAF > ", MAF_threshold," in Chr:", i, " N="),
                                                     paste0("Number of SNPs in", Training_name, "after isolating MAF > ", MAF_threshold," in Chr:", i, " N="),
                                                     paste0("Number of SNPS in ", Training_name, "before flipping Chr:", i," N="),
                                                     paste0("Number of SNPs in ", Validation_name, "before flipping Chr:", i," N="),
                                                     paste0("Combined Number of SNPs BEFORE flipping between ",Validation_name," and ",Training_name," Chr:", i," N="),
                                                     paste0(Validation_name,"_", Training_name, " Chr:", i," remove A-T, C-G Step one: N="),
                                                     paste0(Validation_name,"_", Training_name, " Chr:", i," remove A-T, C-G Step two: N="),
                                                     paste0("Combined all A-T,C-G remove and other mismatches Chr", i," N=")),
                                            N = c(abd[[1]], abd[[2]], abd[[3]], abd[[4]], abd[[5]], abd[[6]], abd[[7]], abd[[8]]))
}
if (INFO_decision == "TRUE"){
  informative.data.frame[[i]] <- data.frame(info = c(paste0("Number of SNPs in", Training_name, "before isolating INFO > ", INFO_threshold," in Chr:", i, " N="),
                                                     paste0("Number of SNPs in", Training_name, "after isolating INFO > ", INFO_threshold," in Chr:", i, " N="),
                                                     paste0("Number of SNPS in ", Training_name, "before flipping Chr:", i," N="),
                                                     paste0("Number of SNPs in ", Validation_name, "before flipping Chr:", i," N="),
                                                     paste0("Combined Number of SNPs BEFORE flipping between ",Validation_name," and ",Training_name," Chr:", i," N="),
                                                     paste0(Validation_name,"_", Training_name, " Chr:", i," remove A-T, C-G Step one: N="),
                                                     paste0(Validation_name,"_", Training_name, " Chr:", i," remove A-T, C-G Step two: N="),
                                                     paste0("Combined all A-T,C-G remove and other mismatches Chr", i," N=")),
                                            N = c(abd[[1]], abd[[2]], abd[[3]], abd[[4]], abd[[5]], abd[[6]], abd[[7]], abd[[8]]))
}
}

for (i in 1:length(chromosomes_to_parse)) {
  filename <- paste0(Validation_name,"_",Training_name,"_dataframe_showing_number_of_SNPs_throughout_PRS_analysis.csv")
  write.table (informative.data.frame[[i]], file = filename, col.names = T, append = T, row.names = F, sep = ",")
}

# data_frames_without_original_row_length <- lapply(informative.data.frame, function(x) {which(x[1,2] == "NA")})
# data_frames_without_original_row_length <- which (data_frames_without_original_row_length == 1)


