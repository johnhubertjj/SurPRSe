### Start Timer
ptm <- proc.time()

#Print R version
getRversion()    

# devtools::install_github('hadley/ggplot2')

library(ggplot2)
library(gridExtra)
library(data.table)

# required packages -- will only install if they are not already installed
list.of.packages <- c("plyr", "stringr", "dplyr", "tidyr", "reshape2", "ggplot2", "scales", "data.table", "plotly","tools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# loads the required packages
lapply(list.of.packages, require, character.only = TRUE)

## set wd
setwd(".")

########################################
# adding in arguments from BASH script #
########################################
args <- commandArgs(trailingOnly = T)
getwd()
print(args)

# Specify the different input tables #
Training_name <- args[3]
Validation_name <- args [4]
Pathway_directory <- args[5]
Pathway_file_name <- args[6] # the name of the file to be accessed (must be in stationary directory)

# This will be a problem, write in a script that defines what the arguments are (probably have to be linked with the arguments script however...stop un-needed analysis)
significance_thresholds <- as.numeric(args[c(7:length(args))])
print(significance_thresholds)



#Training_name <- "CLOZUK_PGC2noclo"
#Validation_name <- "ALSPAC"
#Pathway_directory <- paste0(Training_name,"_",Validation_name,"_output/Pathways/")
#Pathway_file_name <- "/Users/johnhubert/Dropbox/Stationary_data/Pocklington2015_134sets_LoFi_2sets_morphology_notmorphology_deduplicated.txt"
#significance_thresholds <- c(0.05,0.5)

Pathways_used <- file_path_sans_ext(basename(Pathway_file_name))

##### read in pathway sets and standardise column names
pathway_sets <- fread(Pathway_file_name)
setnames(pathway_sets, c("Pathway", "Gene"))

#extract names of each pathway
factorise_column_1<- as.factor(pathway_sets$Pathway)
pathway_names <- levels(factorise_column_1)
number_of_pathways_to_analyse <- length(pathway_names)

final_scoring_output_file <- NULL
order_of_output_for_sig_thresh <- NULL

for (i in significance_thresholds){
scoring_output_file <- paste0(Pathway_directory,pathway_names,"/scoring_",Training_name,"_",Validation_name,"_pathway_",pathway_names,"_with_",i,".profile")
final_scoring_output_file <- c(final_scoring_output_file,scoring_output_file)
order_of_output_for_sig_thresh <- c(order_of_output_for_sig_thresh, paste0(pathway_names,"_",rep(i,length(scoring_output_file))))
}

my_data <- lapply(final_scoring_output_file, read.table, header=TRUE) 
names(my_data) <- str_replace(final_scoring_output_file, pattern = ".profile", replacement = "")


# iterate through significance thresholds 
for (i in 1:length(final_scoring_output_file)) {
  my_data[[i]] <- my_data[[i]][,c(1,2,6)]
  colnames(my_data[[i]]) <- c("FID", "IID", paste("SCORE_", order_of_output_for_sig_thresh[i], sep=""))
}

all_prs <- join_all(my_data, by=c("FID", "IID"), type='left')

write.table(all_prs,file = paste0("FINAL_PATHWAY_RESULTS_PRS_PROFILES",Training_name,"_",Validation_name, "_", Pathways_used, ".txt"),quote = F,row.names = F)
