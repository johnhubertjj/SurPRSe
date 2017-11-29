# Producing significance thresholds arguments flat file # 
# Attempt to produce into same format as plink? # 

#Generate random gene sets from all/brain-expressed genes, where
#probability of selection ~ gene_length, snp_n, indep_snp_n, snp_density, indep_snp_density (didn't use chr)
###################################################################################################

### Start Timer
ptm <- proc.time()

### Read in packages
library (data.table)

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)

# Specify the different input tables #
Training_name <- args[3]
Validation_name <- args[4]
Significance_thresholds <- as.numeric(args[c(5:length(args))])

# Create the plink -q-score-range file for use in plink 
Significance_thresholds_lower_bounds <- unlist(fread(file = paste0(Training_name,"_", Validation_name,"_significance_thresholds_lower_bounds_to_analyse_arguments_file_tmp.txt")))
rownames <- paste0("S",seq(1:length(Significance_thresholds)))
Plink_sig_thresh_file <- data.frame(Significance_thresholds_lower_bounds,Significance_thresholds,row.names = rownames)

# Write out the chromosomes to a text file for easier reading into R in future rather than a bash argument # 
write.table(Plink_sig_thresh_file, file = paste0(Training_name,"_", Validation_name,"_plink_significance_thresholds_arguments_file_tmp.txt"),row.names = T, col.names = F, quote = F)
testing_1 <- fread(paste0(Training_name,"_", Validation_name,"_plink_significance_thresholds_arguments_file_tmp.txt"))      
print(testing_1)

#End Timer
proc.time() - ptm
