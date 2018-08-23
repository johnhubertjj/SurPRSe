#!/bin/env R
#################################################################
##################### PROPER SCRIPT START #######################
#################################################################

### Start Timer
ptm <- proc.time()

library(magrittr)


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
Target_name <- args[1]
Reference_name <- args[2]
Target_dataset <- args[3]
Reference_dataset <- args[4]
Chromosomes <- args[5]

#Reference_dataset <- "/media/johnhubert/PHD DATA/CLOZUK2_noPGC2.assoc.dosage"
#Target_dataset <- "/media/johnhubert/PHD DATA/daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt"

#########################################
# Preparing unique treatment resistance # 
#########################################

Reference_dataset <-"~/Documents/LD_score_COGSandNOCOGS/CLOZUK_PGC2_noCOGS1.tbl"
Target_dataset <- "~/Documents/LD_score_COGSandNOCOGS/daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt"

merged_table <- merge(Reference_table,Target_table,by = "SNP")

Reference_table2 <- merged_table[,.(SNP,BP,CHR,A1.x,A2.x,MAF,OR.x,SE.x,P.x)]
Target_table2 <- merged_table[,.(SNP,BP,CHR,A1.y,A2.y,OR.y,SE.y,P.y)]

setnames(Reference_table2,old = c("A1.x","A2.x","OR.x","SE.x","P.x"), new = c("A1","A2","OR","SE","P"))
setnames(Target_table2, old = c("A1.y","A2.y","OR.y","SE.y","P.y"), new = c("A1","A2","OR", "SE", "P"))

setkey(Reference_table2, CHR)
setkey(Target_table2, CHR)


################################
# Reading in reference dataset #
################################
#hapmap_reffreq_PHASE2_3_CEU <- fread("~/Dropbox/bip1.scz1.ruderfer2014 (1)/hapmap_reffreq_PHASE2_3_CEU.txt", header=T, quote="\"")
thousandgenomes <- fread("~/Documents/1000genomes_european/g1000_European_populations_QC_NO_bad_LD.frq", header=T)
Reference_table <- fread(Reference_dataset)
Target_table <- fread(Target_dataset)

#BP_CLEAN <- merge(BPSCZ.bp.only.results ,hapmap_reffreq_PHASE2_3_CEU,by.x=2,by.y=1)


# drop iINFO out of bounds and refmaf out of bound snps:
#BP_CLEAN_out <- BP_CLEAN[BP_CLEAN[,6] > .9 & BP_CLEAN[,6] < 1.1 & BP_CLEAN[,15] < 0.95 & BP_CLEAN[,15] > 0.05,]


#write.table(BP_CLEAN_out[,1:10], "BIP_CLEAN_FOR_LD_SCORE.txt",row.names=F,quote=F)



for (chromosome.number in 1:22){
  #summ_stat_target<- fread(paste0("CLOZUK2_noPGC2_chr",chromosome.number,".txt"))
  summ_stat_target <- fread(paste0("CLOZUK_PGC2_noCOGS1_table",chromosome.number,".txt"))
  summ_stat_target2 <- copy(summ_stat_target)
  #Training_table_1 <- fread(paste0("./RSupdate/CLOZUK_chr",chromosome.number,".KAVIAR.update"))
  summ_stat_reference <- fread(paste0("CLOZUK_PGC2noclo_table",chromosome.number,".txt"))
  summ_stat_reference_2 <- copy(summ_stat_reference)
  
  LD_score_reference <- fread("w_hm3.snplist")
  LD_score_reference_2 <- read.table(paste0("eur_w_ld_chr/",chromosome.number,".l2.ldscore.gz"),stringsAsFactors = F,header = T)
  
  #colnames(Training_table_1) <- c("SNP","RS_SNP")
  #summ_stat_target_RS <- merge(Training_table_1, summ_stat_target, by="SNP", all = F)
  #summ_stat_reference[, RS_SNP := SNP ]
  #CHR <- args[8] 
  #SNP <- args[9]
  #BP <- args[10]
  #A1 <- args[11]
  #A2 <- args[12]
  #OR <- args[13]
  #BETA <- args[14]
  #number_of_frequency_columns <- args[16]
  
  
  # match SNPs with common identifiers 
  #combined_table_1 <- merge(summ_stat_target, summ_stat_reference, by = c("CHR","BP","A1","A2"), all = F)
  
  # alter the allele positions
  #altered_table_reference_1 <- setnames(summ_stat_reference, old = c("A1","A2"), new = c("A2","A1"))
  
  #merge again
  #combined_table_2 <- merge(summ_stat_target,altered_table_reference_1, by = c("CHR","BP","A1","A2"), all = F)
  #combined_table_3 <- merge(combined_table_1,combined_table_2, by = c("CHR","BP"), all = F)
  
  #final_table <- rbind(combined_table_1, combined_table_2)
  
  #chartr(old = c("A","T","C","G"), new = c("T", "G", "C", "A")) 
  
  
  
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
  log.to.odds <- function(imported.data.table,name) {
    imported.dt.col.names <- colnames(imported.data.table)
    if (any("OR" == imported.dt.col.names) == F) {
      cat("No Odds Ratio included in", deparse(substitute(imported.data.table)))
    }else{
      testing <- deparse(substitute(imported.data.table))
      assign(paste0(name,".BETA"), log(imported.data.table$OR), envir = e)
      imported.data.table$BETA <- log(imported.data.table$OR)
      assign(testing, imported.data.table, envir = e)
    }
  }
  
  Beta.to.odds <- function(imported.data.table, name) {
    imported.dt.col.names <- colnames(imported.data.table)
    if (any("BETA" == imported.dt.col.names) == F) {
      cat("No BETA coefficient included in", deparse(substitute(imported.data.table)))
    }else{
      testing <- deparse(substitute(imported.data.table))
      assign(paste0(name,".OR"), exp(imported.data.table$BETA), envir = e)
      imported.data.table$OR <- exp(imported.data.table$BETA)
      assign(testing, imported.data.table, envir = e)
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
  PGC.duplicated.removed <- which(duplicated(summ_stat_target$SNP))
  PGC.duplicated.removed.rev <- which(duplicated(summ_stat_target$SNP, fromLast = T))
  
  ### One-line duplication check - Training ###
  length(summ_stat_target$SNP)
  if (length(PGC.duplicated.removed) >= 1){
    PGC_duplicate_SNPS <- summ_stat_target$SNP[(duplicated(summ_stat_target$SNP) | duplicated(summ_stat_target$SNP, fromLast = TRUE)) ]
    summ_stat_target <- summ_stat_target[!(duplicated(summ_stat_target$SNP) | duplicated(summ_stat_target$SNP, fromLast = TRUE)) ]
  }
  length(summ_stat_target$SNP)
  
  cat("Number of SNPS in ", Target_name," after duplication check Chr:",chromosome.number, "N=" ,nrow(summ_stat_target))
  
  
  cat("Number of SNPS in ", Reference_name,"_chromosome_", chromosome.number,": N=" ,nrow(summ_stat_reference))
  
  ### Remove Duplicated SNPs in Training here ###
  PGC.duplicated.removed <- which(duplicated(summ_stat_reference$SNP))
  PGC.duplicated.removed.rev <- which(duplicated(summ_stat_reference$SNP, fromLast = T))
  
  ### One-line duplication check - Training ###
  length(summ_stat_reference$SNP)
  if (length(PGC.duplicated.removed) >= 1){
    PGC_duplicate_SNPS <- summ_stat_reference$RS_SNP[(duplicated(summ_stat_reference$RS_SNP) | duplicated(summ_stat_reference$SNP, fromLast = TRUE)) ]
    summ_stat_reference <- summ_stat_reference[!(duplicated(summ_stat_reference$RS_SNP) | duplicated(summ_stat_reference$SNP, fromLast = TRUE)) ]
  }
  length(summ_stat_reference$SNP)
  
  cat("Number of SNPS in ",Reference_name," after duplication check Chr:",chromosome.number, "N=" ,nrow(summ_stat_reference))
  
  Beta.to.odds(summ_stat_target, Target_name)
  log.to.odds(summ_stat_target, Target_name)
  Beta.to.odds(summ_stat_reference, Reference_name)
  log.to.odds(summ_stat_reference, Reference_name)
  
  
  # Simple tracking of each row for future checks
  e$summ_stat_target$place <- c(1:length(summ_stat_target$CHR))
  e$summ_stat_reference$place <- c(1:length(summ_stat_reference$CHR))
  
  ### Merging based on BP position ###
  Combined_summ_stat <- merge(x = e$summ_stat_reference, y = e$summ_stat_target, by = c('BP','CHR'), all = F, sort = F)
  Combined_summ_stat <- Combined_summ_stat[,RS_SNP:= SNP.y]
  #setnames(LD_score_reference, old = c("SNP","A1","A2"), new = c("RS_SNP","A1_LD","A2_LD"))
  #setnames(LD_score_reference_2, old = c("CHR","SNP","BP"), new = c("CHR_LD","RS_SNP","BP_LD"))
  #Combined_summ_stat <- merge(x = Combined_summ_stat, y = LD_score_reference, by = "RS_SNP")
  #Combined_summ_stat <- merge(x = Combined_summ_stat, y = LD_score_reference_2, by = "RS_SNP")
  
 # Change all alleles to uppercase: saving time method
  # below it states:
  ## assign columns you want to change to object cols
  ## take all rows of data.table; apply the toupper function to Subset of Data columns named in cols object
  cols <- c('A1.x','A2.x','A1.y','A2.y',"A1_LD","A2_LD")
  Combined_summ_stat[, (cols) := lapply(.SD, toupper), .SDcols = cols ]
  
  ### Checking and limiting to SNPs with only one allele for each A1 and A2
  Checking_length_of_alleles(Combined_summ_stat)
  
  new.target.datatable <- Combined_summ_stat[,.(RS_SNP,BP,CHR,A1.y,A2.y,OR.y,SE.y,P.y)]
  setnames(new.target.datatable, old = c("RS_SNP","A1.y","A2.y","OR.y","SE.y","P.y"), new = c("SNP","A1","A2","OR", "SE", "P"))
  
  new.reference.datatable <- Combined_summ_stat[,.(RS_SNP,BP,CHR,A1.x,A2.x,OR.x,SE.x,P.x)]
  setnames(new.reference.datatable, old = c("RS_SNP","A1.x","A2.x","OR.x","SE.x","P.x"), new = c("SNP","A1","A2","OR", "SE", "P"))
  
  new.target.table <- paste0(Target_name,"_table_changed_no_flipped_alleles", chromosome.number,"_new.txt")
  new.reference.table <-  paste0(Reference_name,"_table_changed_no_flipped_alleles", chromosome.number,"_new.txt")
  
  # Write update file for plink
  write.table(new.target.datatable, file = new.target.table, quote = F, row.names = F, col.names = T)
  
  # Write update file for plink
  write.table(new.reference.datatable, file = new.reference.table, quote = F, row.names = F, col.names = T)
}
  
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
  
  # Change the values of the OR to suit Allele changes #
  if (length(a) >= 1) {
    # Combined_summ_stat$BETA[a] <- (-Combined_summ_stat$BETA[a])
    PGC.OR <- Combined_summ_stat$OR.x[a]
    #PGC.BETA <- Combined_summ_stat$BETA[a]
    change.odds(PGC.OR)
    #change.beta(PGC.BETA)
    Combined_summ_stat$OR.x[a]<- e$PGC.NEW.OR
    #Combined_summ_stat$BETA[a] <- e$PGC.NEW.BETA
    rm(PGC.NEW.OR, envir = e)
    rm(PGC.OR)
    #rm(PGC.NEW.BETA, envir = e)
    #rm(PGC.BETA)
    Combined_summ_stat$A1.x[a] <- as.character(Combined_summ_stat$A1.y[a])
    Combined_summ_stat$A2.x[a] <- as.character(Combined_summ_stat$A2.y[a])
  }
  
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
  
  # Same OR change #
  if (length(a) >= 1){
    # Combined_summ_stat$BETA[a] <- (-Combined_summ_stat$BETA[a])
    PGC.OR <- Combined_summ_stat$OR.x[a]
    #PGC.BETA <- Combined_summ_stat$BETA[a]
    change.odds(PGC.OR)
    #change.beta(PGC.BETA)
    Combined_summ_stat$OR.x[a]<- e$PGC.NEW.OR
    #Combined_summ_stat$BETA[a] <- e$PGC.NEW.BETA
    rm(PGC.NEW.OR, envir = e)
    rm(PGC.OR)
    #rm(PGC.NEW.BETA, envir = e)
    #rm(PGC.BETA)
    
    Combined_summ_stat$A1.x[a] <- as.character(Combined_summ_stat$A1.y[a])
    Combined_summ_stat$A2.x[a] <- as.character(Combined_summ_stat$A2.y[a])
  }
  
  a <- which(Combined_summ_stat$A1.y != Combined_summ_stat$A1.x | Combined_summ_stat$A2.y!=Combined_summ_stat$A2.x)
  if (length(a)>0) Combined_summ_stat <- Combined_summ_stat[-a,]
  cat("Chr:", chromosome.number,", remove A-T, C-G and other mismatches: N=", nrow(Combined_summ_stat))
  
  ### Find  out where in the original table the SNP identifiers are
  PGC.final.index <- Combined_summ_stat$place.x
  PGC.A1.allele.changes <- Combined_summ_stat$A1.x
  PGC.A2.allele.changes <- Combined_summ_stat$A2.x
  PGC.OR.changes <- Combined_summ_stat$OR.x
  #PGC.BETA.changes <- Combined_summ_stat$BETA
  CLOZUK.final.index <- Combined_summ_stat$place.x
  
  
  ### Alter original files to have the right allele switching and the right BETA coefficient #
  summ_stat_reference <- summ_stat_reference[PGC.final.index, A1:= PGC.A1.allele.changes]
  summ_stat_reference <- summ_stat_reference[PGC.final.index,A2:= PGC.A2.allele.changes]
  summ_stat_reference <- summ_stat_reference[PGC.final.index,OR:= PGC.OR.changes]
  # summ_stat_reference <- summ_stat_reference[PGC.final.index,BETA:= PGC.BETA.changes]
  oldnames_training <- colnames(summ_stat_reference)
  
  
  ### checking everything is above board ###
  Checking_allele_swapping(summ_stat_reference, summ_stat_reference_2, which.is.combined = "NONE")
  Checking_allele_swapping(summ_stat_target, summ_stat_target2, which.is.combined = "NONE")
  
  if (length(e$summ_stat_targetflipped.alleles) != 0) {
    warning("alleles have been flipped on CLOZUK")
  }
  checking.which.OR.do.not.equal.each.other <- which(summ_stat_reference$OR != summ_stat_reference_2$OR)
  
  problematic.SNPs <- match(e$summ_stat_referenceflipped.alleles, checking.which.OR.do.not.equal.each.other)
  problematic.SNPs <- which(is.na(problematic.SNPs))
  check1_index <- e$summ_stat_referenceflipped.alleles[problematic.SNPs]
  check1 <- summ_stat_reference[check1_index]
  check2 <- summ_stat_reference[check1_index]
  
  if (all(check1$OR == 1) & all(check2$OR == 1)){
    cat("You have",length(check1$CHR),"flipped SNPs with OR = 1")
  } else {
    warning("There is an uneven amount of flipped SNPs")
  }
  
  ### Replacing PGC values###
  PGC.integers.to.change <- Combined_summ_stat[,.I[grep(paste0("^",chromosome.number,":\\d+:\\w+:\\w+"),SNP.x, perl = T, invert = T)]]
  PGC.divisive.names <- Combined_summ_stat[PGC.integers.to.change]
  Combined_summ_stat <- Combined_summ_stat[,SNP.x := paste0(CHR,":",BP)]
  
  PGC.integers.to.change <- Combined_summ_stat[,.I[grep(paste0("^",chromosome.number,":\\d+:\\w+:\\w+"),SNP.y, perl = T, invert = T)]]
  PGC.divisive.names <- Combined_summ_stat[PGC.integers.to.change]
  Combined_summ_stat <- Combined_summ_stat[,SNP.y := paste0(CHR,":",BP)]
  

  new.target.datatable <- Combined_summ_stat[,.(RS_SNP,BP,CHR,A1.y,A2.y,OR.y,SE.y,P.y)]
  setnames(new.target.datatable, old = c("RS_SNP","A1.y","A2.y","OR.y","SE.y","P.y"), new = c("SNP","A1","A2","OR", "SE", "P"))
  
  new.reference.datatable <- Combined_summ_stat[,.(RS_SNP,BP,CHR,A1.x,A2.x,OR.x,SE.x,P.x)]
  setnames(new.reference.datatable, old = c("RS_SNP","A1.x","A2.x","OR.x","SE.x","P.x"), new = c("SNP","A1","A2","OR", "SE", "P"))
  
  new.LD.SNPlist <- Combined_summ_stat[,.(RS_SNP,A1_LD,A2_LD)]
  setnames(new.LD.SNPlist, old = c("RS_SNP","A1_LD","A2_LD"), new = c("SNP","A1","A2"))
  
  new.LD.score_reference <- Combined_summ_stat[,.(CHR_LD,RS_SNP,BP_LD,CM,MAF,L2)]
  setnames(new.LD.score_reference, old = c("RS_SNP","CHR_LD","BP_LD"), new = c("SNP","CHR","BP"))
  
  new.target.table <- paste0(Target_name,"_table_changed", chromosome.number,"_new.txt")
  new.reference.table <-  paste0(Reference_name,"_table_changed", chromosome.number,"_new.txt")
  new.LD_snp_list_ref <- paste("LD_SNPlist_chromosome_",chromosome.number,"_new.txt")
  new.LD.score_reference.table <- paste0("eur_w_ld_chr/analagous_LD_score_references",chromosome.number,"_new.txt")
  
  # Write target file 
  write.table(new.target.datatable, file = new.target.table, quote = F, row.names = F, col.names = T)
  
  # Write reference file
  write.table(new.reference.datatable, file = new.reference.table, quote = F, row.names = F, col.names = T)
  
  # Write LD SNPlist 
  write.table(new.LD.SNPlist, file = new.LD_snp_list_ref, quote = F, row.names = F, col.names = T)
  
  # Write LD table 
  write.table(new.LD.score_reference, file = new.LD.score_reference.table, quote = F, row.names = F, col.names = T)
}
  
## Read in the latest Summary stats tables after converting to one table
for (i in 1:22){
  assign(paste0("COGS", i), fread(paste0("CLOZUK_PGC2noclo_table_changed_no_flipped_alleles",i,"_new.txt")),envir = .GlobalEnv)
}
l = list()

## print out to one table under a common filename
for (i in 1:22) {
  l[[i]] <- eval(parse(text = paste0("COGS",i)))
}
combined_final_table <- rbindlist(l)

write.table(combined_final_table, file = "CLOZUK_PGC2noclo_full_table_aligned_and_matched_snps.txt",quote = F,row.names = F)

## Read in the latest Summary stats tables after converting to one table
for (i in 1:22){
  assign(paste0("NO_COGS", i), fread(paste0("CLOZUK_PGC2_noCOGS_table_changed_no_flipped_alleles",i,"_new.txt")),envir = .GlobalEnv)
}
l = list()

## print out to one table under a common filename
for (i in 1:22) {
  l[[i]] <- eval(parse(text = paste0("NO_COGS",i)))
}
combined_final_table <- rbindlist(l)

write.table(combined_final_table, file = "CLOZUK_PGC2_noCOGS_full_table_aligned_and_matched_snps.txt",quote = F,row.names = F)



## Read in the latest Summary stats tables after converting to one table
for (i in 1:22){
  assign(paste0("NON_trs_table", i), fread(paste0("Non_TRS_people_aka_PGC2_table",i,"_new.txt")),envir = .GlobalEnv)
}
l = list()

## print out to one table under a common filename
for (i in 1:22) {
  l[[i]] <- eval(parse(text = paste0("NON_trs_table",i)))
}
combined_final_table <- rbindlist(l)

write.table(combined_final_table, file = "PGC_GWAS_metal_aligned_and_matched_snps.txt",quote = F,row.names = F)

## Read in the latest Summary stats tables after converting to one table
for (i in 1:22){
  assign(paste0("TRS_table", i), fread(paste0("TRS_people_aka_CLOZUK_table",i,"_new.txt")),envir = .GlobalEnv)
}
l = list()

## print out to one table under a common filename
for (i in 1:22) {
  l[[i]] <- eval(parse(text = paste0("NON_trs_table",i)))
}
combined_final_table <- rbindlist(l)

write.table(combined_final_table, file = "CLOZUK.1kGp1_table.txt",quote = F,row.names = F)
###########################################
###### GWIS #######
###########################################

h2TRS <- 0.3549 #(0.0183)
h2non_TRS <-  0.4942 #(0.0255)
gcov <-  0.4188 #(0.0216)

Non_TRS <- fread("PGC_GWAS_metal_aligned_and_matched_snps.txt")
TRS <-  fread("CLOZUK.1kGp1_table.txt")
thousandgenomes <- fread("~/Documents/1000genomes_european/g1000_European_populations_QC_NO_bad_LD.frq", header=T)

BPSCZ.bp.only.results <- read.delim("~/Dropbox/bip1.scz1.ruderfer2014 (1)/BPSCZ.bp-only.results.txt",stringsAsFactors =F)
BPSCZ.scz.only.results <- read.delim("~/Dropbox/bip1.scz1.ruderfer2014 (1)/BPSCZ.scz-only.results.txt",stringsAsFactors =F)


#### Rotate to obtain Unique BIP and Unique SCZ, first merge the files and allign the reference allele:
primary <- merge(TRS,Non_TRS,by="SNP")

#primary  <- merge(BPSCZ.bp.only.results,BPSCZ.scz.only.results,by=2)

primary[primary$A1.x != primary$A1.y,7] <- exp(log(primary[primary$A1.x != primary$A1.y,7]) * -1)

primary[primary$A1.x != primary$A1.y,4] <- primary[primary$A1.x != primary$A1.y,13]

primary[primary$A2.x != primary$A2.y,5] <- primary[primary$A2.x != primary$A2.y,14]
View(primary)


# merge with reference allel frequency:
primary2 <- merge(primary,thousandgenomes,by="SNP")

# flip ref allele freq to align with the ref allel in BIP and SCZ
rows_with_MAF_change <- which(primary2$A1.x != primary2$A1) 
primary2$MAF[rows_with_MAF_change] <- 1 - primary2$MAF[rows_with_MAF_change]


# compute allel frequencies to ensure proper accounting for the AF in computing the new beta:
p <- primary2$MAF
base <- 2*p * (1-p) + (1-p)^2
het <- 2*p * (1-p) / base
hom <- (1-p)^2/base




################################# BELOW FIRST BIP - SCZ #################################
# TRS - NON_TRS


primary2 <- primary2[,1:24]

# compute the variables c as described in the methosds, essential for transformation.
c <- (h2TRS - gcov) / (h2TRS + gcov)


primary2 <- cbind(primary2, het * ( (1+ c) *log(primary2$OR.x) - (1 - c) * log(primary2$OR.y))/2 + hom * ( (1+ c) * 2*log(primary2$OR.x)- (1 - c) * 2*log(primary2$OR.y))/4 )

primary2 <- primary2[,1:25]

# USe the delta method to get the SE for the beta:
require(msm)
se <- matrix(NA,nrow=nrow(primary2),ncol=1)
for(i in 1:nrow(primary2)){
  c1 <- 1+c 
  c2 <- 1-c 
  homi <- hom[i]
  heti <- het[i]
  se[i,] <-  deltamethod(~(heti*(c1) *log(x2) - (c2) *log(x1))/2 +  homi*((c1) *2*log(x2) - (c2) * 2*log(x1))/4 ,mean=c(primary2$OR.y[i],primary2$OR.x[i]),cov=diag(c(primary2$SE.y[i],primary2$SE.x[i])) %*%matrix(c(1,0,0,1),2,2) %*% diag(c(primary2$SE.y[i],primary2$SE.x[i])) )
} 

primary2 <- cbind(primary2,se)
primary2 <- cbind(primary2,pchisq((primary2[,25]/primary2[,26])^2,1,lower.tail=F))


primary2_out <- primary2[primary2[,6] > .9 & primary2[,6] < 1.1 & primary2[,24] < 0.95 & primary2[,24] > 0.05,]


BIP_min_SCZ_out <- primary2_out[,c(1,2,3,4,5,24,25,26,27)]

names(BIP_min_SCZ_out) <- c("SNP","CHR", "BP","A1","A2","EAF","Beta","SE","P" )

write.table(BIP_min_SCZ_out,"BIP_min_SCZ_out_LDSCORE.txt",row.names=F,quote=F)


######################################## REPEAT FOR SCZ - BIP #############################################

# drop all BIP steps

primary2 <- primary2[,1:24]



# compute the variables c as described in the methosds, essential for transformation.
c <- (h2BIP - gcov) / (h2BIP + gcov)


primary2 <- cbind(primary2, het * ( (1+c) *log(primary2$OR.y) -  (1-c) * log(primary2$OR.x))/2 + hom * ( (1+ c) * 2*log(primary2$OR.y)- (1 - c) * 2*log(primary2$OR.x))/4 )

primary2 <- primary2[,1:25]

# Use the delta method to get the SE for the beta:
require(msm)
se <- matrix(NA,nrow=nrow(primary2),ncol=1)
for(i in 1:nrow(primary2)){
  c1 <- 1+c 
  c2 <- 1-c 
  homi <- hom[i]
  heti <- het[i]
  se[i,] <-  deltamethod(~(heti*(c1) *log(x1) - (c2) *log(x2))/2 +  homi*((c1) *2*log(x1) - (c2) * 2*log(x2))/4 ,mean=c(primary2$OR.y[i],primary2$OR.x[i]),cov=diag(c(primary2$SE.y[i],primary2$SE.x[i])) %*%matrix(c(1,0,0,1),2,2) %*% diag(c(primary2$SE.y[i],primary2$SE.x[i])) )
} 

primary2 <- cbind(primary2,se)
primary2 <- cbind(primary2,pchisq((primary2[,25]/primary2[,26])^2,1,lower.tail=F))


primary2_out <- primary2[primary2[,6] > .9 & primary2[,6] < 1.1 & primary2[,24] < 0.95 & primary2[,24] > 0.05,]


SCZ_min_BIP_out <- primary2_out[,c(1,2,3,4,5,24,25,26,27)]

names(SCZ_min_BIP_out) <- c("SNP","CHR", "BP","A1","A2","EAF","Beta","SE","P" )

write.table(SCZ_min_BIP_out,"SCZ_min_BIP_out_LDSCORE.txt",row.names=F,quote=F)














