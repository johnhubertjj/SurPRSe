#!/bin/env R
#################################################################
##################### PROPER SCRIPT START #######################
#################################################################

### Start Timer
ptm <- proc.time()

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)

# MAGMA making a SNP loc file

library(data.table)
library(tools)

# Function for GMT creator

writeGMT <- function #Create a gmt (gene matrix transposed) file
### Createss a gmt (gene matrix transposed) file such as those
### provided by mSigDB or geneSigDB, from an R list object.
### Function by Levi Waldron.
(object,
 ### R list object that will be converted to GMT file.  Each element
 ### should contain a vector of gene names, and the names of the
 ### elements will used for the gene set names
 fname
 ### Output file name for .gmt file
){
  if (class(object) != "list") stop("object should be of class 'list'")
  if(file.exists(fname)) unlink(fname)
  for (iElement in 1:length(object)){
    write.table(t(c(make.names(rep(names(object)[iElement],2)),object[[iElement]])),
                sep="\t",quote=FALSE,
                file=fname,append=TRUE,col.names=FALSE,row.names=FALSE)
  }
  ### Called for the effect of writing a .gmt file
}

Training_name <- args[3]
Validation_name <- args[4]
Gene_set_file <- args[5]
input_file_Path <- paste0("combined_", Training_name,"_", Validation_name, "_table_with_CHR.POS_identifiers_flipped_alleles_no_duplicates_MAF_INFO.txt")

input_file <- fread(input_file_Path)
SNP_loc <- input_file[,.(SNP,CHR,BP)]

training_name_literal <- file_path_sans_ext(input_file_Path)

write.table(SNP_loc, file = paste0(training_name_literal,"snploc.txt"),quote = F,row.names = F,col.names = T)


# Reading in Andrews annotations


annot <- fread(Gene_set_file)

basename_annot <- basename(Gene_set_file)
annot_name <- file_path_sans_ext(basename_annot)

list_of_gmt <- list()
rownames <- unique(annot$V1)
for (i in rownames){
  genes <- annot[V1 == i,V2]
  list_item <- c(genes)
  list_of_gmt[[i]] <- list_item
}

writeGMT(list_of_gmt,fname = paste0(annot_name,".gmt"))
