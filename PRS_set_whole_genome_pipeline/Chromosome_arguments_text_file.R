# Producing Chromosomes argument flat file # 

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
Training_name <- args[3]
Validation_name <- args[4]
Chromosomes_to_analyse <- as.numeric(args[c(5:length(args))])

# Write out the chromosomes to a text file for easier reading into R in future rather than a bash argument # 
write(Chromosomes_to_analyse, file = paste0(Training_name, "_", Validation_name,"_chromosomes_to_analyse_arguments_file_tmp.txt"),ncolumns = 1)
testing_1 <- scan(file = paste0(Training_name, "_", Validation_name,"_chromosomes_to_analyse_arguments_file_tmp.txt"),sep="\n")      
print(testing_1)

#End Timer
proc.time() - ptm
