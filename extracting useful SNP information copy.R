## Load packages ## 
library(stringr)
library(data.table)
library(plyr)
### gain the number of files that match the chromosomes used ###
setwd("/Users/johnhubert/Documents/PR54/PGC_CLOZUK_PRS/extrainfo/Rout_files/")
current_files <- system("ls",intern = T)
current_files <- c(current_files,1,"ahroghouhguo[1]")
current_files_index <- grep("CLOZUK_PGC_COMBINE_chr\\d+.Rout",current_files,perl = T)
current_files <- current_files[current_files_index]

### For loop (change to apply to batch script "PBS_NUM_ARRAY") ####
file_collection <- list()
for (i in 1:length(current_files)) {
  file_collection[[i]] <- read.delim(paste0("CLOZUK_PGC_COMBINE_chr",i,".Rout"))
} 

### Find the rows where cat function appears ###
informative.data.frame <- list()

for (i in 1:length(current_files)) {
length_of_PGC_CLOZUK_before_integer <- grep("Read\\s\\d+\\srows",file_collection[[i]][,1],perl = T)
length_of_combined_datasets_integer <- grep("Chr\\:\\s\\d+",file_collection[[i]][,1],perl = T)
new_vector_one <- as.character(file_collection[[i]][length_of_PGC_CLOZUK_before_integer,1])
new_vector_two <- as.character(file_collection[[i]][length_of_combined_datasets_integer,1])
abc <-str_extract_all(new_vector_one,"(?<=Read\\s)\\d+(?=\\srows)")
abd <-str_extract_all(new_vector_two,"\\d+(?=\\>)")
if (length(unlist(abc)) != 2) {
  informative.data.frame[[i]] <- data.frame(info = c(paste0("PGC chr",i," N="),
                                                     paste0("CLOZUK chr",i," N="),
                                                     paste0("Combined chr",i," N="),
                                                     paste0("Combined first AT-CG remove chr",i," N="),
                                                     paste0("Combined all A-T,C-G remove chr",i," N="),
                                                     paste0("Combined all A-T,C-G remove and other mismatches chr",i," N=")),
                                            N = c("NA","NA",abd[[1]],abd[[2]],abd[[3]],abd[[4]]))
}else{
informative.data.frame[[i]] <- data.frame(info = c(paste0("PGC chr",i," N="),
                                           paste0("CLOZUK chr",i," N="),
                                           paste0("Combined chr",i," N="),
                                           paste0("Combined first AT-CG remove chr",i," N="),
                                           paste0("Combined all A-T,C-G remove chr",i," N="),
                                           paste0("Combined all A-T,C-G remove and other mismatches chr",i," N=")),
                                     N = c(abc[[1]],abc[[2]],abd[[1]],abd[[2]],abd[[3]],abd[[4]]))
}
}


data_frames_without_original_row_length <- lapply(informative.data.frame, function(x) {which(x[1,2] == "NA")})
data_frames_without_original_row_length <- which (data_frames_without_original_row_length == 1)

for (i in data_frames_without_original_row_length) {
  setwd("/Users/johnhubert/Dropbox/PGC/")
  PGC.test.data.frame <- fread(paste0("gzip -dc PGC_table",i,".txt.gz"))
  informative.data.frame[[i]][,2] <- as.character(informative.data.frame[[i]][,2])
  informative.data.frame[[i]][1,2] <- nrow(PGC.test.data.frame)
  
  setwd("/Users/johnhubert/Documents/BGE/")
  untar(paste0("CLOZUK_GWAS_BGE_chr",i,".tar.gz"),files = paste0("CLOZUK_GWAS_BGE_chr",i,".bim"))
  CLOZUK.data <- fread(paste0("CLOZUK_GWAS_BGE_chr",i,".bim"))
  informative.data.frame[[i]][2,2] <- nrow(CLOZUK.data)
  rm(PGC.test.data.frame)
  rm(CLOZUK.data)
}

total_PGC_SNPs <- lapply(informative.data.frame, function (x) {
  as.character(x[1,2])
})
total_PGC_SNPs <- sum(as.numeric(total_PGC_SNPs))

total_CLOZUK_SNPs <- lapply(informative.data.frame, function (x) {
  as.character(x[2,2])
})
total_CLOZUK_SNPs <- sum(as.numeric(total_CLOZUK_SNPs))

total_combined_SNPs <- lapply(informative.data.frame, function (x) {
  as.character(x[3,2])
})
total_combined_SNPs <- sum(as.numeric(total_combined_SNPs))

total_number_of_SNPs_after_comparison <- lapply(informative.data.frame, function (x) {
  as.character(x[6,2])
})
total_number_of_SNPs_after_comparison <- sum(as.numeric(total_number_of_SNPs_after_comparison))

percent.SNP.remaining <- (total_number_of_SNPs_after_comparison / total_combined_SNPs)*100
percent.SNP.lost <- ((total_combined_SNPs - total_number_of_SNPs_after_comparison) / total_combined_SNPs) * 100

informative.data.frame[[23]] <- data.frame(info = c("Total PGC SNPs =",
                                                    "Total CLOZUK SNPs =",
                                                    "Total Combined SNPs =",
                                                    "Total Number of SNPs after removal =",
                                                    "Percentage SNPs lost =",
                                                    "Percentage SNPs remain ="),
                                            N =   c(as.character(total_PGC_SNPs),
                                                    as.character(total_CLOZUK_SNPs),
                                                    as.character(total_combined_SNPs),
                                                    as.character(total_number_of_SNPs_after_comparison),
                                                    as.character(percent.SNP.lost),
                                                    as.character(percent.SNP.remaining)))

setwd(setwd("/Users/johnhubert/Documents/PR54/PGC_CLOZUK_PRS/extrainfo/Rout_files/"))

save.image(file="CLOZUK_PGC_combination_info.RData",compress=T)

informative.list.to.df <- ldply(informative.data.frame, function(t) t$toDataFrame())

write.table (informative.data.frame[[1]],file = "informative_data_frame.csv", col.names = T, row.names = F, sep = ",") 

for (i in 2:23) {
write.table (informative.data.frame[[i]],file = "informative_data_frame.csv", col.names = T, append = T, row.names = F, sep = ",")
}

