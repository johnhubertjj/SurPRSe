#!/bin/env R
#################################################################
##################### PROPER SCRIPT START #######################
#################################################################

### Start Timer
ptm <- proc.time()



## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite("Biostrings")
library(Biostrings)
library(data.table)

#######################################
# adding in arguments from BASH script#
#######################################
args <- commandArgs(trailingOnly = T)
print(args)


# specify the different input tables #
Target_name <- args[3]
Reference_name <- args[4]
Target_dataset <- args[5]
Reference_dataset <- args[6]
Chromosomes <- args[7]

Reference_dataset <- "/media/johnhubert/PHD DATA/CLOZUK2_noPGC2.assoc.dosage"
Target_dataset <- "/media/johnhubert/PHD DATA/daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt"

summ_stat_target<- fread("CLOZUK2_noPGC2_chr1.txt")
Training_table_1 <- fread("./RSupdate/CLOZUK_chr1.KAVIAR.update")
summ_stat_reference <- fread("PGC2noCLOZUK_chr1.txt")

colnames(Training_table_1) <- c("SNP","RS_SNP")
summ_stat_target_RS <- merge(Training_table_1, summ_stat_target, by="SNP", all = F)
summ_stat_reference[, RS_SNP := SNP ]

#CHR <- args[8] 
#SNP <- args[9]
#BP <- args[10]
#A1 <- args[11]
#A2 <- args[12]
#OR <- args[13]
#BETA <- args[14]
#number_of_frequency_columns <- args[16]


# match SNPs with common identifiers 
combined_table_1 <- merge(summ_stat_target, summ_stat_reference, by = c("CHR","BP","A1","A2"), all = F)

# alter the allele positions
altered_table_reference_1 <- setnames(summ_stat_reference, old = c("A1","A2"), new = c("A2","A1"))

#merge again
combined_table_2 <- merge(summ_stat_target,altered_table_reference_1, by = c("CHR","BP","A1","A2"), all = F)
combined_table_3 <- merge(combined_table_1,combined_table_2, by = c("CHR","BP"), all = F)

final_table <- rbind(combined_table_1, combined_table_2)

chartr(old = c("A","T","C","G"), new = c("T", "G", "C", "A")) 
       
       
       
##################################################################
#### CHECK TO SEE IF TWO COLUMNS HAVE ALLELES SWAPPED FUNCTION ###
##################################################################
       
       Checking_allele_swapping <- function(alteredtable1,table2,which.is.combined = c("NONE","PGC","CLOZUK")){
         if (which.is.combined == "PGC") {
           Allele.x <- alteredtable1[,c("A1.x","A2.x"),with = F]
           Allele.y <- table2[,c("A1","A2"), with = F]
         }
         if(which.is.combined == "CLOZUK") {
           Allele.x <- alteredtable1[,c("A1.y","A2.y"),with = F]
           Allele.y <- table2[,c("A1","A2"), with = F]
         }
         if(which.is.combined == "NONE"){
           Allele.x <- alteredtable1[,c("A1","A2"),with = F]
           Allele.y <- table2[,c("A1","A2"), with = F]
         }
         
         Allele.x <- as.data.frame(Allele.x)
         Allele.y <- as.data.frame(Allele.y)
         iterations.to.remain <- NULL
         
         for (i in 1:length(Allele.x[,1])){
           if (Allele.x[i,1] != Allele.y[i,1]) {
             iterations.to.remain <- c(iterations.to.remain,i)
           }
         }
         argument.name <-deparse(substitute(alteredtable1))
         assign(paste0(argument.name,"flipped.alleles"), iterations.to.remain, envir = e)
       }
       
       Checking_length_of_alleles <- function(input_table1){
         # Removes all multiple allele counts (aka indels and deletions) to clean up the data.
         input_table1 <- input_table1[!(nchar(A1.x) > 1 | nchar(A2.x) > 1)]
         input_table1 <- input_table1[!(nchar(A1.y) > 1 | nchar(A2.y) > 1)]
         
         # Checks for any Alleles which are not standard 1:1
         integers_to_check <- which(input_table1$A1.x != "A" & input_table1$A1.x != "C" & input_table1$A1.x !=  "T" & input_table1$A1.x !=  "G")
         integers_to_check2 <- which(input_table1$A1.y != "A" & input_table1$A1.y != "C" & input_table1$A1.y !=  "T" & input_table1$A1.y !=  "G")
         integers_to_check3 <- which(input_table1$A2.x != "A" & input_table1$A2.x != "C" & input_table1$A2.x !=  "T" & input_table1$A2.x !=  "G")
         integers_to_check4 <- which(input_table1$A2.y != "A" & input_table1$A2.y != "C" & input_table1$A2.y !=  "T" & input_table1$A2.y !=  "G")
         
         # Would probably work better in a loop but I got too lazy
         if(length(integers_to_check) != 0){
           input_table1 <- input_table1[!integers_to_check]
         }
         
         if(length(integers_to_check2) != 0){
           input_table1 <- input_table1[!integers_to_check2]
         }
         
         if(length(integers_to_check3) != 0){
           input_table1 <- input_table1[!integers_to_check3]
         }
         
         if(length(integers_to_check4) != 0){
           input_table1 <- input_table1[!integers_to_check4]
         }
         assign("Combined_summ_stat", input_table1, envir = .GlobalEnv)
       }

       ###############################
       # ODDS RATIO TO BETA FUNCTION #
       ###############################
       # add new environment#
       e <- new.env()
       
       # BETA and OR change functions ###
       log.to.odds <- function(imported.data.table) {
         imported.dt.col.names <- colnames(imported.data.table)
         if (any("OR" == imported.dt.col.names) == F) {
           cat("No Odds Ratio included in", deparse(substitute(imported.data.table)))
         }else{
           assign("PGC.BETA", log(imported.data.table$OR), envir = e)
         }
       }
       
       Beta.to.odds <- function(imported.data.table) {
         imported.dt.col.names <- colnames(imported.data.table)
         if (any("BETA" == imported.dt.col.names) == F) {
           cat("No BETA coefficient included in", deparse(substitute(imported.data.table)))
         }else{
           assign("PGC.OR", exp(imported.data.table$BETA), envir = e)
         }
       }
       
       change.odds <- function (odds.ratios) {
         PGC.NEW.OR <- 1 / odds.ratios
         assign("PGC.NEW.OR", PGC.NEW.OR, envir = e)
       }
       
       change.beta <- function (beta.coefficients) {
         PGC.NEW.BETA <- -(beta.coefficients)
         assign("PGC.NEW.BETA", PGC.NEW.BETA, envir = e)
       }
####################################################################
       cat("Number of SNPS in ",Target_name,"_chromosome_", chromosome.number,": N=" , nrow(summ_stat_target))
       
       ### Remove Duplicated SNPs in Training here ###
       PGC.duplicated.removed <- which(duplicated(summ_stat_target$RS_SNP))
       PGC.duplicated.removed.rev <- which(duplicated(summ_stat_target$RS_SNP, fromLast = T))
       
       ### One-line duplication check - Training ###
       length(summ_stat_target$RS_SNP)
       if (length(PGC.duplicated.removed) >= 1){
         PGC_duplicate_SNPS <- summ_stat_target$SNP[(duplicated(summ_stat_target$RS_SNP) | duplicated(summ_stat_target$RS_SNP, fromLast = TRUE)) ]
         summ_stat_target <- summ_stat_target[!(duplicated(summ_stat_target$RS_SNP) | duplicated(summ_stat_target$RS_SNP, fromLast = TRUE)) ]
       }
       length(summ_stat_target$RS_SNP)
       
       cat("Number of SNPS in ", Target_name," after duplication check Chr:",chromosome.number, "N=" ,nrow(summ_stat_target))
       
       
       cat("Number of SNPS in ", Reference_name,"_chromosome_", chromosome.number,": N=" ,nrow(summ_stat_reference))
       
       ### Remove Duplicated SNPs in Training here ###
       PGC.duplicated.removed <- which(duplicated(summ_stat_reference$RS_SNP))
       PGC.duplicated.removed.rev <- which(duplicated(summ_stat_reference$RS_SNP, fromLast = T))
       
       ### One-line duplication check - Training ###
       length(summ_stat_reference$SNP)
       if (length(PGC.duplicated.removed) >= 1){
         PGC_duplicate_SNPS <- summ_stat_reference$RS_SNP[(duplicated(summ_stat_reference$RS_SNP) | duplicated(summ_stat_reference$RS_SNP, fromLast = TRUE)) ]
         summ_stat_reference <- summ_stat_reference[!(duplicated(summ_stat_reference$RS_SNP) | duplicated(summ_stat_reference$RS_SNP, fromLast = TRUE)) ]
       }
       length(summ_stat_reference$RS_SNP)
       
       cat("Number of SNPS in ",Reference_name," after duplication check Chr:",chromosome.number, "N=" ,nrow(summ_stat_reference))
       
       
       
       
       # Simple tracking of each row for future checks
       summ_stat_target$place <- c(1:length(summ_stat_target$CHR))
       summ_stat_reference$place <- c(1:length(summ_stat_reference$CHR))
       
       ### Merging based on BP position ###
       Combined_summ_stat <- merge(x = summ_stat_target_RS, y = summ_stat_reference, by = c('BP','CHR'), all = F, sort = F)
       
       ### Checking and limiting to SNPs with only one allele for each A1 and A2
       Checking_length_of_alleles(Combined_summ_stat)
       
       # Change all alleles to uppercase: saving time method
       # below it states:
       ## assign columns you want to change to object cols
       ## take all rows of data.table; apply the toupper function to Subset of Data columns named in cols object
       cols <- c('A1.x','A2.x','A1.y','A2.y')
       Combined_summ_stat[, (cols) := lapply(.SD, toupper), .SDcols = cols ]
       
       cat("Number of SNPs BEFORE flipping between CLOZUK and PGC Chr:", chromosome.number , "N=", nrow(Combined_summ_stat))
       
      # if (nrow(Combined_summ_stat) > nrow(PGC.alternative) | nrow(Combined_summ_stat) > nrow(CLOZUK.alternative)){
      #   warning("combined dataset size is larger than one/both of the input datasets, check for duplicates")
      # }
       
       a <- which(Combined_summ_stat$A1.x == "C" & Combined_summ_stat$A2.x == "G"); if (length(a)>0) {Combined_summ_stat<-Combined_summ_stat[-a]}
       a <- which(Combined_summ_stat$A1.x == "G" & Combined_summ_stat$A2.x == "C"); if (length(a)>0) {Combined_summ_stat<-Combined_summ_stat[-a]}
       a <- which(Combined_summ_stat$A1.x == "A" & Combined_summ_stat$A2.x == "T"); if (length(a)>0) {Combined_summ_stat<-Combined_summ_stat[-a]}
       a <- which(Combined_summ_stat$A1.x == "T" & Combined_summ_stat$A2.x == "A"); if (length(a)>0) {Combined_summ_stat<-Combined_summ_stat[-a]}
       
       cat(Target_name,"_", Reference_name, "Chr:",chromosome.number,"remove A-T, C-G Step one: N=" ,nrow(Combined_summ_stat))
       
       a <- which(Combined_summ_stat$A1.y=="C" & Combined_summ_stat$A2.y=="G"); if (length(a)>0) {Combined_summ_stat<-Combined_summ_stat[-a]}
       a <- which(Combined_summ_stat$A1.y=="G" & Combined_summ_stat$A2.y=="C"); if (length(a)>0) {Combined_summ_stat<-Combined_summ_stat[-a]}
       a <- which(Combined_summ_stat$A1.y=="A" & Combined_summ_stat$A2.y=="T"); if (length(a)>0) {Combined_summ_stat<-Combined_summ_stat[-a]}
       a <- which(Combined_summ_stat$A1.y=="T" & Combined_summ_stat$A2.y=="A"); if (length(a)>0) {Combined_summ_stat<-Combined_summ_stat[-a]}
       
       
       cat("CLOZUK_PGC Chr:", chromosome.number, "remove A-T, C-G Step two: N=", nrow(Combined_summ_stat))
       
       ## FLIPPING ##
       # A1.x/A2.x are Training set alleles
       # A1.y/A2.y are Genotype set alleles
       
       a <- which (Combined_summ_stat$A1.y == Combined_summ_stat$A2.x & Combined_summ_stat$A2.y == Combined_summ_stat$A1.x)
       
       # Find alleles that are swapped in a more complicated manner (eg: with A -> G)
       a <- which(Combined_summ_stat$A1.y == Combined_summ_stat$A1.x & Combined_summ_stat$A2.y == Combined_summ_stat$A2.x)
       d <- seq(1:nrow(Combined_summ_stat));d <- d[-a]
       
       # Specify which alleles in dataset 1 have A1 = A and A2 = G
       test1 <- which (Combined_summ_stat$A1.x[d] == "A" & Combined_summ_stat$A2.x[d] == "G")
       test1a <- which (Combined_summ_stat$A1.x[d] == "G" & Combined_summ_stat$A2.x[d] == "A")
       test1 <- unique(c(test1, test1a))
       
       if (length(test1) > 0){
         for (i in 1:length(test1)){
           if(Combined_summ_stat$A1.y[d[test1[i]]] == "T" & Combined_summ_stat$A2.y[d[test1[i]]] == "C"){
             if(Combined_summ_stat$A1.x[d[test1[i]]] == "A"){
               Combined_summ_stat$A1.x[d[test1[i]]] <- "T"
               Combined_summ_stat$A2.x[d[test1[i]]] <- "C"
             }
           }else if (Combined_summ_stat$A1.y[d[test1[i]]] == "C" & Combined_summ_stat$A2.y[d[test1[i]]] == "T"){
             if(Combined_summ_stat$A1.x[d[test1[i]]] == "G") {
               Combined_summ_stat$A1.x[d[test1[i]]] <- "C"
               Combined_summ_stat$A2.x[d[test1[i]]] <- "T"
             }
           }else{
             Combined_summ_stat[-d[test1[i]]]
           }
         }
       }
       
       
       test2 <- which(Combined_summ_stat$A1.x[d] == "T" & Combined_summ_stat$A2.x[d] == "C")
       test2a <- which(Combined_summ_stat$A1.x[d] == "C" & Combined_summ_stat$A2.x[d] == "T")
       test2 <- unique(c(test2, test2a))
       
       if (length(test2) > 0){
         for (i in 1:length(test2)){
           if(Combined_summ_stat$A1.y[d[test2[i]]] != "A" & Combined_summ_stat$A2.y[d[test2[i]]] != "G" | Combined_summ_stat$A1.y[d[test2[i]]] != "G" & Combined_summ_stat$A2.y[d[test2[i]]] != "A"){
             Combined_summ_stat[-d[test2[i]]]
           }else{
             if(Combined_summ_stat$A1.x[d[test2[i]]] == "T"){
               Combined_summ_stat$A1.x[d[test2[i]]] <- "A"
               Combined_summ_stat$A2.x[d[test2[i]]] <- "G"
             }
             if(Combined_summ_stat$A1.x[d[test2[i]]] == "C") {
               Combined_summ_stat$A1.x[d[test2[i]]] <- "G"
               Combined_summ_stat$A2.x[d[test2[i]]] <- "A"
             }
           }
         }
       }
       
       a <- which(Combined_summ_stat$A1.y==Combined_summ_stat$A2.x & Combined_summ_stat$A2.y==Combined_summ_stat$A1.x)
       
       a <- which(Combined_summ_stat$A1.y != Combined_summ_stat$A1.x | Combined_summ_stat$A2.y!=Combined_summ_stat$A2.x)
       if (length(a)>0) Combined_summ_stat <- Combined_summ_stat[-a,]
       cat("Chr:", chromosome.number,", remove A-T, C-G and other mismatches: N=", nrow(Combined_summ_stat))
       
       ### Find  out where in the original table the SNP identifiers are
       PGC.final.index <- Combined_summ_stat$place.x
       PGC.A1.allele.changes <- Combined_summ_stat$A1.x
       PGC.A2.allele.changes <- Combined_summ_stat$A2.x
       PGC.OR.changes <- Combined_summ_stat$OR
       PGC.BETA.changes <- Combined_summ_stat$BETA
       CLOZUK.final.index <- Combined_summ_stat$place.x
       
       
       ### Alter original files to have the right allele switching and the right BETA coefficient #
       PGC.data.frame <- PGC.data.frame[PGC.final.index,A1:= PGC.A1.allele.changes]
       PGC.data.frame <- PGC.data.frame[PGC.final.index,A2:= PGC.A2.allele.changes]
       PGC.data.frame <- PGC.data.frame[PGC.final.index,OR:= PGC.OR.changes]
       PGC.data.frame <- PGC.data.frame[PGC.final.index,BETA:= PGC.BETA.changes]
       oldnames_training <- colnames(PGC.data.frame)
       
  