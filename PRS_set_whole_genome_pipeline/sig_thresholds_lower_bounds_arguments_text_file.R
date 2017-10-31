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

# specify the different input tables #
training_set_name <- args[3]
validation_set_name <- args[4]
Significance_thresholds_lower_bounds <- as.numeric(args[c(5:length(args))])

# Write out the chromosomes to a text file for easier reading into R in future rather than a bash argument # 
write(Significance_thresholds_lower_bounds, file = paste0(training_set_name,"_", validation_set_name,"_significance_thresholds_lower_bounds_to_analyse_arguments_file_tmp.txt"), ncolumns = 1)
testing_1 <- fread(file = paste0(training_set_name,"_", validation_set_name,"_significance_thresholds_lower_bounds_to_analyse_arguments_file_tmp.txt"))      
print(testing_1)

#End Timer
proc.time() - ptm
