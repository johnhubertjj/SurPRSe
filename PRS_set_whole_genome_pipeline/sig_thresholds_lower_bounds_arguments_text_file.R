# Producing significance thresholds arguments flat file # 
# Attempt to produce into same format as plink? # 

#Generate random gene sets from all/brain-expressed genes, where
#probability of selection ~ gene_length, snp_n, indep_snp_n, snp_density, indep_snp_density (didn't use chr)
###################################################################################################

### Start Timer
ptm <- proc.time()

### Read in required packages
#list.of.packages <- c("data.table")

### Find packages which are not installed on this version of R 
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#Install these packates
#if(length(new.packages)) install.packages(new.packages)

# Find old packages
#old_packages <- old.packages()
#old_packages_use <- as.data.frame(old_packages,stringsAsFactors = F)

# Find out if any of the currently installed packages are old and install the newer, stable versions
#if(any(list.of.packages %in% old_packages_use$Package)){
#  df_of_packages <- data.frame(list.of.packages,stringsAsFactors = F)
#  colnames(df_of_packages) <- "Package"
#  packages_to_update_table <- merge(df_of_packages,old_packages, by = "Package", all = F)
#  Update_packages <- packages_to_update_table$Package
#  install.packages(Update_packages)
#}

# loads the required packages
#lapply(list.of.packages, require, character.only = TRUE)

library("data.table")

tmp.install.packages <- function(pack, dependencies=TRUE, ...) {
  path <- tempdir()
  ## Add 'path' to .libPaths, and be sure that it is not
  ## at the first position, otherwise any other package during
  ## this session would be installed into 'path'
  firstpath <- .libPaths()[1]
  .libPaths(c(firstpath, path))
  install.packages(pack, dependencies=dependencies, lib=path, ...)
}
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
testing_1 <- fread(paste0(training_set_name,"_", validation_set_name,"_significance_thresholds_lower_bounds_to_analyse_arguments_file_tmp.txt"))      
print(testing_1)

#End Timer
proc.time() - ptm
