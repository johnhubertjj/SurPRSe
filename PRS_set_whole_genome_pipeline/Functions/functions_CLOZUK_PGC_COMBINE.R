# CLOZUK_PGC_COMBINE_functions #

###############################
# ODDS RATIO TO BETA FUNCTION #
###############################
# add new environment#
e <- new.env()

# initialise summary information table
df <- tibble(No_info1 = NA, No_info2 = NA, No_info3 = NA, No_info4 = NA)
assign("SNP.information.table", df, envir = e)

### BETA and OR change functions ###
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
  assign("combined.CLOZUK.PGC", input_table1, envir = .GlobalEnv)
}

###############################################
#### Record number of SNPs in each section ####
###############################################

Add_information <- function(df, string_of_information,e, stage, chromosome){
  
  if(all(grepl("No_info",colnames(e$SNP.information.table)) ) == TRUE){
    current_tibble <- tibble(x = string_of_information, y = nrow(df), stage = stage, chromosome = chromosome)
    assign("SNP.information.table", current_tibble, envir = e)
    
  }else{
    current_tibble <- add_row(e$SNP.information.table, x = string_of_information, y = nrow(df), stage = stage, chromosome = chromosome)
    assign("SNP.information.table", current_tibble, envir = e)
    }
}

