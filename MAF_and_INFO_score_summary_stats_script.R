### Start Timer
ptm <- proc.time()

### Library
library(data.table)

###############################
# Checking location for serial or batch analysis #
###############################

System_information <- Sys.info()
whereami <- System_information['user']

if (whereami == "johnhubert") {
  chromosome.number <- 22
  a <- "PGC_table22"
  
} else if (whereami == 'JJ') {
  chromosome.number <- 22
  a <- "PGC_table22"
  
} else if (whereami == "c1020109") {
  
  # Preparing to run in Job array
  AI <- Sys.getenv("PBS_ARRAY_INDEX")
  chromosome.number <- as.numeric(AI)
  
} else {
  stop("current environment NOT at work/home or on servers, please add to script above to specify where you are and what type of analysis you want to do")
}

#######################################
# adding in arguments from BASH script#
#######################################
args=(commandArgs(trailingOnly = T))

if (length(args) == 0) {
  print ("No arguments supplied, assumed to be local analysis")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

# specify the different input tables #
Training_datatable <- paste0("./output/",a,".txt")
Training_datatable_output <- paste0("./output/",a,"_new.txt")
Training_name <- c
MAF_decision <- e
INFO_decision <- f

# Read in Training set
PGC_table <- fread(Training_datatable)
changed_PGC_table <- copy(PGC_table)

if (INFO_decision == "YES") {
# remove all SNPs with INFO > 0.9
changed_PGC_table <- changed_PGC_table[INFO > 0.9]
}

# Do I want to remove SNPs based on MAF
if (MAF_decision == "YES") {
  
# remove all SNPs lower that 0.01
# change so that the column recognises other datasets
FRQ_to_keep_helpful3 <- which(changed_PGC_table$FRQ_A_29415 <= 0.01 & changed_PGC_table$FRQ_U_40101 <= 0.01)
changed_PGC_table <- changed_PGC_table[!FRQ_to_keep_helpful3]

#Find out which columns have one value below MAF 0.01
FRQ_U_helpful2 <- which(changed_PGC_table$FRQ_A_29415 <= 0.01 & changed_PGC_table$FRQ_U_40101 > 0.01)
FRQ_A_helpful1 <- which(changed_PGC_table$FRQ_A_29415 > 0.01 & changed_PGC_table$FRQ_U_40101 <= 0.01)

# add extra column with MAF for each SNP according to the rule, keep the minimum allele frequency (except those below 0.01)
changed_PGC_table[,MAF := pmin(FRQ_A_29415,FRQ_U_40101)]
changed_PGC_table[FRQ_A_helpful1, MAF := changed_PGC_table[FRQ_A_helpful1,FRQ_A_29415]]
changed_PGC_table[FRQ_U_helpful2, MAF := changed_PGC_table[FRQ_U_helpful2,FRQ_U_40101]]

# do a check 
if (all(changed_PGC_table$MAF == changed_PGC_table$FRQ_A_29415 | changed_PGC_table$FRQ_U_40101) == FALSE){
  warning("MAF is not neccessarily correct for this chromosome, check back please")
}
}

# write new table overwriting the previous PGC_new_table
write.table(changed_PGC_table, file = Training_datatable_output, row.names = F, quote = F)
cat("Number of SNPs in", Training_name, "before isolating MAF > 0.01 and INFO > 0.9 =", nrow(PGC_table))
cat("Number of SNPs in", Training_name, "after isolating MAF > 0.01 and INFO > 0.9 =", nrow(changed_PGC_table))

#End Timer
proc.time() - ptm