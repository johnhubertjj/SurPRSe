### Start Timer
ptm <- proc.time()

## New environment
e <- new.env()

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
Training_datatable <- paste0("./", Training_name, "_", Validation_name,"_output/", args[3],"_new.txt")
Validation_datatable_bim_file <- paste0(args[4],".bim")
Training_name <- args[3]
Validation_name <- args[4]
Rout_file_storage <- args[5]
INFO_decision <- args[6]
MAF_decision <- args[7]
MAF_threshold <- args[8]
INFO_threshold <- args[9]
SE_decision <- args[10]
SE_threshold <- args[11]
chromosomes_to_parse <- as.numeric(args[c(12:length(args))])
print(Chromosomes_to_parse)

#test<- function(Training_datatable, Validation_datatable_bim_file, Training_name, Validation_name, Rout_file_storage, INFO_decision, INFO_threshold, MAF_decision, MAF_threshold, SE_decision, SE_threshold, chromosomes_to_parse){
# Deciding the first line of the csv file function
Print_thresholds <- function(INFO_decision, MAF_decision, SE_decision){
  whoishere <- rep(FALSE, 3)
    if (MAF_decision == TRUE){
    whoishere[1]  <- TRUE
    }
  
  if(INFO_decision == TRUE){
    whoishere[2]  <- TRUE
  }

  if(SE_decision == TRUE){
    whoishere[3]  <- TRUE  
  }
  whoishere <- paste(whoishere, collapse = " ")
  assign("whoishere", whoishere, envir = e)
}

# writing the first line of csv file function
Figuring_out_first_two_lines <- function(whoishere, INFO_threshold, MAF_threshold, SE_threshold){
  if(grepl("TRUE TRUE TRUE", whoishere)){
    line <- paste0("MAF > ", MAF_threshold, ", INFO > ", INFO_threshold, "and SE > ", SE_threshold)
  }
  if(grepl("TRUE TRUE FALSE",whoishere)){
    line <- paste0("MAF > ", MAF_threshold, " and INFO > ", INFO_threshold)
  }
  if(grepl("TRUE FALSE FALSE", whoishere)){
    line <- paste0("MAF > ", MAF_threshold)
  }
  if(grepl("FALSE FALSE FALSE", whoishere)){
    line <- "no thresholds"
  }
  if(grepl("FALSE TRUE FALSE", whoishere)){
    line <- paste0("INFO > ", INFO_threshold)
  }
  if(grepl("FALSE TRUE TRUE", whoishere)){
    line <- paste0("INFO > ", INFO_threshold, "and SE > ", SE_threshold)
  }
  if(grepl("FALSE FALSE TRUE", whoishere)){
    line <- paste0("SE > ", SE_threshold)
  }
  if(grepl("TRUE FALSE TRUE", whoishere)){
    line <-  line <- paste0("MAF > ", MAF_threshold, " and SE > ", SE_threshold)
  }
  assign("line", line, envir = e)
}
  

### gain the number of files that match the chromosomes used ###
setwd(Rout_file_storage)
current_files <- system("ls", intern = T)


### For loop for file parsing ####
file_collection <- list()
for (i in chromosomes_to_parse) {
  current_files_index <- grep(paste0("\\d+\\[",i,"\\].raven0.OU"),current_files,perl = T)
  file_collection[[i]] <- read.delim(current_files[current_files_index])
} 

### Find the rows where cat function appears
informative.data.frame <- list()
Warnings_here_terminal <- list()
Errors_here_terminal <- list()


for (i in chromosomes_to_parse) {
  browser()
# initially collect all the interesting information # 
length_of_combined_datasets_integer <- grep("Chr\\:\\s\\d+",file_collection[[i]][,1],perl = T)
all_parse <- grep("N=\\s\\d+", file_collection[[i]][,1],perl = T)

new_vector_two <- as.character(file_collection[[i]][length_of_combined_datasets_integer,1])
new_vector_one <- as.character(file_collection[[i]][all_parse,1])

### Check for warning messages that have been flagged:terminal
warnings <- grep("Warning\\:", file_collection[[i]][,1], perl = T)

if(length(warnings) > 0){
Warnings_here_terminal[[i]] <- as.character(file_collection[[i]][warnings,1])
}

# Catching errors from the terminal
errors <- grep("Error\\:", file_collection[[i]][,1], perl = T)

if(length(errors) > 0){
  Errors_here_terminal[[i]] <- as.character(file_collection[[i]][errors,1])
}

# Extract all strings with a space before a number at the end : should only be a chromosome number for now
abc <- str_extract_all(new_vector_one,"(?<=\\s)\\d+$")
abd <-str_extract_all(new_vector_two,"(?<=\\s)\\d+$")

# Get the right print-out
Print_thresholds(INFO_decision,MAF_decision, SE_decision)
Figuring_out_first_two_lines(e$whoishere, INFO_threshold, MAF_threshold, SE_threshold) 

# Print out the number of SNPs throughout the pipeline

  informative.data.frame[[i]] <- data.frame(info = c(paste0("Number of SNPs in ", Training_name, " before analysis"),
                                                     paste0("Number of SNPs in ", Training_name, " after FIRST duplication check"),
                                                     paste0("Number of SNPs in ", Training_name, " before isolating ", e$line, " in Chr:", i, " N="),
                                                     paste0("Number of SNPs in ", Training_name, " after isolating ", e$line," in Chr:", i, " N="),
                                                     paste0("Number of SNPS in ", Training_name, " before flipping Chr:", i," N="),
                                                     paste0("Number of SNPs in ", Training_name, " after SECOND duplication check"),
                                                     paste0("Number of SNPs in ", Validation_name, " before flipping Chr:", i," N="),
                                                     paste0("Number of SNPs in ", Validation_name, " after FIRST duplication check"),
                                                     paste0("Combined Number of SNPs BEFORE flipping between ", Validation_name," and ", Training_name," Chr:", i," N="),
                                                     paste0(Validation_name,"_", Training_name, " Chr:", i," remove A-T, C-G Step one: N="),
                                                     paste0(Validation_name,"_", Training_name, " Chr:", i," remove A-T, C-G Step two: N="),
                                                     paste0("Combined all A-T,C-G remove and other mismatches Chr", i," N=")),
                                            N = c(abc[[1]], abc[[2]], abc[[3]], abc[[4]], abc[[5]], abc[[6]], abc[[7]], abc[[8]], abc[[9]], abc[[10]], abc[[11]], abc[[12]]),stringsAsFactors = F)
if(length(warnings) > 0){
  for (warn in 1:length(warnings)){
    informative.data.frame[[i]] <- rbind(informative.data.frame[[i]], c(Warnings_here_terminal[[i]][warn], ""))  
  }
} 

if(length(errors) > 0){
  for (err in 1:length(errors)){
    informative.data.frame[[i]] <- rbind(informative.data.frame[[i]], c(Errors_here_terminal[[i]][err], ""))  
  }
}
  
}
#}

Final_data_frame <- data.frame(info = c(paste0("Number of SNPs in ", Training_name, " before analysis"),
                                        paste0("Number of SNPs in ", Training_name, " after FIRST duplication check"),
                                        paste0("Number of SNPs in ", Training_name, " before isolating ", e$line, ": N="),
                                        paste0("Number of SNPs in ", Training_name, " after isolating ", e$line,": N="),
                                        paste0("Number of SNPS in ", Training_name, " before flipping: N="),
                                        paste0("Number of SNPs in ", Training_name, " after SECOND duplication check"),
                                        paste0("Number of SNPs in ", Validation_name, " before flipping: N="),
                                        paste0("Number of SNPs in ", Validation_name, " after FIRST duplication check"),
                                        paste0("Combined Number of SNPs BEFORE flipping between ", Validation_name," and ", Training_name,": N="),
                                        paste0(Validation_name,"_", Training_name, " remove A-T, C-G Step one: N="),
                                        paste0(Validation_name,"_", Training_name, " remove A-T, C-G Step two: N="),
                                        paste0("Combined all A-T,C-G remove and other mismatches: N=")),
                               N = rep(0,12))

current_sum <- rep(0,length(chromosomes_to_parse))

for (length in 1:12){
  for (i in chromosomes_to_parse) {
    current_sum[i] <- informative.data.frame[[i]][length,2]
  }
current_sum <- as.numeric(current_sum)
final_sum <- sum(current_sum)
Final_data_frame[length,2] <- final_sum
Final_data_frame[length,1] <- informative.data.frame[[1]][length,1]
}

combined_data_frame <- rbindlist(informative.data.frame)


filename_whole_genome <- paste0(Validation_name,"_",Training_name,"_dataframe_showing_number_of_SNPs_throughout_PRS_analysis_whole_genome.csv")
filename_by_chromosome <- paste0(Validation_name,"_",Training_name,"_dataframe_showing_number_of_SNPs_throughout_PRS_analysis_by_chromosome.csv")
  
write.table (combined_data_frame, file = filename_by_chromosome, col.names = T, append = T, row.names = F, sep = ",")
write.table (Final_data_frame, file = filename_whole_genome, col.names = T, append = T, row.names = F, sep = ",")
  

# data_frames_without_original_row_length <- lapply(informative.data.frame, function(x) {which(x[1,2] == "NA")})
# data_frames_without_original_row_length <- which (data_frames_without_original_row_length == 1)


