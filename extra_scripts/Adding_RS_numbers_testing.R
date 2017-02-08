###EXPERIMENTATION###



### Checking to see if two columns have alleles swapped ###

Checking_allele_swapping <- function(Merged_table){
  Allele.x <- Merged_table[,c("A1.x","A2.x"),with = F]
  Allele.y <- Merged_table[,c("A1.y","A2.y"),with = F]
  
  Allele.x <- as.data.frame(Allele.x)
  Allele.y <- as.data.frame(Allele.y)
  iterations.to.remain <- NULL
  
  for (i in 1:length(Allele.x[,1])){
    if ((((Allele.x[i,1] == Allele.y[i,1]) | (Allele.x[i,1] == Allele.y[i,2])) & ((Allele.x[i,2] == Allele.y[i,1]) | (Allele.x[i,2] == Allele.y[i,2])))){
      iterations.to.remain <- c(iterations.to.remain,i)
    }
  }
  Merged_table2 <<- Merged_table[iterations.to.remain]  
}






# eval(parse) is used to evaluate an R expression from a string, so because you are looping through objects, eval(parse) makes it possible...

setnames(CLOZUK.data, c("CHR","SNP","BETA","BP","A1","A2"))

# Match the rs Numbers with that of CLOZUK
# Will take into account no matches and keep the previous information
common.alleles.index <- match(CLOZUK.data$BP, chr22$BP)
CLOZUK.data.2 <- CLOZUK.data$SNP
indexed_Na_positions <- which(is.na(common.alleles.index))  
rs.numbers.with.problems  <- grep("rs.*:\\d+:.*:.*", CLOZUK.data$SNP, perl = T)

## WRITE A DIFFERENCE BETWEEN SNP Values from CLOZUK.data and based on BP positions...
## Alleles are different at certain locations, so may be excluded anyway but still need to check...
CLOZUK.data.3 <- CLOZUK.data.2
CLOZUK.data.2 <- chr22[common.alleles.index,BP]
CLOZUK.data.2 <- paste0("rs", CLOZUK.data.2)
CLOZUK.data.2[indexed_Na_positions] <- CLOZUK.data.3[indexed_Na_positions]

rs.numbers.with.problems  <- grep("rs.*:\\d+:.*:.*", CLOZUK.data.2, perl = T)
splitting_rs_number_strings <- CLOZUK.data.2[rs.numbers.with.problems]
multiple_splitting<-strsplit(splitting_rs_number_strings, ":")
Multiple_splitting_2 <- sapply(X = multiple_splitting, FUN = "[", 1)
CLOZUK.data.2[rs.numbers.with.problems] <- Multiple_splitting_2
CLOZUK.data$RS <- CLOZUK.data.2

################################################################################
# checking which rs numbers are duplicated, all are right next to each other#

index_duplicated_data.CLOZUK <- which(duplicated(CLOZUK.data$BP) == T)
index_duplicated_data1.CLOZUK <- which(duplicated(CLOZUK.data$BP,fromLast = T) == T)
index_duplicated_data.CLOZUK <- c(index_duplicated_data.CLOZUK,index_duplicated_data1.CLOZUK)
CLOZUK_data_duplicated.CLOZUK <- CLOZUK.data[index_duplicated_data.CLOZUK]

index_duplicated_data.PGC <- which(duplicated(PGC_test_data_frame$BP) == T)
index_duplicated_data1.PGC <- which(duplicated(PGC_test_data_frame$BP,fromLast = T) == T)
index_duplicated_data.PGC <- c(index_duplicated_data.PGC,index_duplicated_data1.PGC)
PGC_data_duplicated.PGC <- PGC_test_data_frame[index_duplicated_data.PGC]

PGC_duplicated_merge <- merge (CLOZUK_data_duplicated.CLOZUK[,c("CHR","SNP","BP","A1","A2"),with = F ],PGC_data_duplicated.PGC[,c("CHR","SNP","BP","A1","A2"), with = F],by = c('BP'), all = F)
PGC_dupicated_merge2 <- apply (PGC_duplicated_merge,1,function(x) {})
################################################################################

# how many identifiers are left which are unknown?
sum(!grepl ("rs",CLOZUK.data.2))  
length(unique(CLOZUK.data.2))


#Matching rs values

## checking overlap between both CLOZUK and schizophrenia data
PGC_integer_positions <- match(PGC_test_data_frame$SNP,CLOZUK.data$RS)
CLOZUK_integer_positions <- match(CLOZUK.data$RS, PGC_test_data_frame$SNP)

indexed_Na_positions.CLOZUK <- which(is.na(CLOZUK_integer_positions))  
indexed_Na_positions.PGC <- which(is.na(PGC_integer_positions))

PGC.final.positions <- which(!is.na(PGC_integer_positions))
CLOZUK.final.positions <- which(!is.na(CLOZUK_integer_positions))

Combined.PGC.CLOZUK.PGC.results <- PGC_test_data_frame[PGC.final.positions,]
Combined.PGC.CLOZUK.CLOZUK.results <- CLOZUK.data[CLOZUK.final.positions,]

test2 <- match(Combined.PGC.CLOZUK.CLOZUK.results$RS,Combined.PGC.CLOZUK.PGC.results$SNP)
test3 <- which(!is.na(test2))

PGC.merge.df <- PGC_test_data_frame[,1:5, with = F]
CLOZUK.merge.df <- CLOZUK.data[,c(1:2,4:6),with = F]
length(unique(test4$SNP))

test4 <- merge (PGC.merge.df,CLOZUK.merge.df,by = c('CHR','SNP','BP'), all = F)

a <- CLOZUK.merge.df[c(524009,524010),]
b <- PGC.merge.df[148774]

###replacing rs values PGC###
common.alleles.index.PGC <- match(PGC_test_data_frame$BP, chr22$BP)
PGC.data.2 <- PGC_test_data_frame$SNP
indexed_Na_positions.PGC.2 <- which(is.na(common.alleles.index.PGC))  
PGC.data.3 <- PGC.data.2
PGC.data.2 <- chr22[common.alleles.index.PGC,RS]
PGC.data.2 <- paste0("rs", PGC.data.2)
PGC.data.2[indexed_Na_positions.PGC.2] <- PGC.data.3[indexed_Na_positions.PGC.2]
PGC_test_data_frame$SNP <- PGC.data.2

PGC_test_data_frame2 <- PGC_test_data_frame
PGC_test_data_frame3 <- rep(0,nrow(PGC_test_data_frame2[,'CHR']))


### Testing to see if you can use fread to read in files quicker than usual...####
read.entry.table <- function(file, entry) {
  browser()
  lines <- readLines(file)
  
  table.entry <- lines == entry
  if (sum(table.entry) != 1) stop(paste(entry, "not found"))
  
  empty.lines <- which(lines == "")
  empty.lines <- c(empty.lines, length(lines) + 1L)
  
  table.start <- which(table.entry) + 1L
  table.end   <- empty.lines[which(empty.lines > table.start)[1]] - 1L
  
  return(read.table(textConnection(lines[seq(from = table.start,
                                             to   = table.end)]),
                    header = TRUE))
}


# create sample dataset
set.seed(1)
m   <- matrix(rnorm(1e5),ncol=10)
csv <- data.frame(x=1:1e4,m)
write.csv(csv,"test.csv")
# s: rows we want to read
s <- c(1:50,53, 65,77,90,100:200,350:500, 5000:6000)
# v: logical, T means read this row (equivalent to your read_vec)
v <- (1:1e4 %in% s)

seq  <- rle(v)
idx  <- c(0, cumsum(seq$lengths))[which(seq$values)] + 1
# indx: start = starting row of sequence, length = length of sequence (compare to s)
indx <- data.frame(start=idx, length=seq$length[which(seq$values)])

library(data.table)
result <- do.call(rbind,apply(indx,1, function(x) return(fread("test.csv",nrows=x[2],skip=x[1]))))

#### Separating out the PGC_data
#### 
PGC_data <- fread(paste(fpath,"PGC/daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt",sep = ""))

for (i in 22:1) {
  assign(paste0("PGC_chr",i), PGC_data[,.I[grep(i,CHR)]])
  assign(paste0("PGC.chr.table",i), PGC_data[eval(parse(text = paste0("PGC_chr",i)))])
  z <- gzfile(paste0(fpath,"PGC/PGC_table",i,".txt.gz"))
  write.table(eval(parse(text = paste0("PGC.chr.table",i))), file = z ,row.names = F,quote = F)
  a <- c(paste0("PGC_chr",i),paste0("PGC.chr.table",i))
  rm(list = a)
}
##### 
#data.table manipulation
set.seed(1)
DT <- data.table(station=sample.int(n=3, size=1e6, replace=TRUE), 
                 wind=rgamma(n=1e6, shape=1.5, rate=1/10),
                 other=rnorm(n=1.6),
                 key="station")
idx = DT [, .I[which(c(NA, diff(wind)) > 35)], by=station][, V1]
DT[idx, wind := NA_real_]

######
flipping_function <- function() {
  browser()
  a <- which (combined.CLOZUK.PGC$A1.y == combined.CLOZUK.PGC$A2.x & combined.CLOZUK.PGC$A2.y == combined.CLOZUK.PGC$A1.x)
  # combined.CLOZUK.PGC$BP[a] <- (-combined.CLOZUK.PGC$BP[a])
  combined.CLOZUK.PGC$A1.x[a] <- as.character(combined.CLOZUK.PGC$A1.y[a])
  combined.CLOZUK.PGC$A2.x[a] <- as.character(combined.CLOZUK.PGC$A2.y[a])
  
  a <- which(combined.CLOZUK.PGC$A1.y == combined.CLOZUK.PGC$A1.x & combined.CLOZUK.PGC$A2.y == combined.CLOZUK.PGC$A2.x)
  d <- seq(1:nrow(combined.CLOZUK.PGC));d <- d[-a]
  
  b1 <- which(combined.CLOZUK.PGC$A1.x[d] == "A")
  b2 <- which(combined.CLOZUK.PGC$A1.x[d] == "C")
  b3 <- which(combined.CLOZUK.PGC$A1.x[d] == "G")
  b4 <- which(combined.CLOZUK.PGC$A1.x[d] == "T")
  
  combined.CLOZUK.PGC$A1.x[d[b1]] <- "T"
  combined.CLOZUK.PGC$A1.x[d[b2]] <- "G"
  combined.CLOZUK.PGC$A1.x[d[b3]] <- "C"
  combined.CLOZUK.PGC$A1.x[d[b4]] <- "A"
  
  b1 <- which(combined.CLOZUK.PGC$A2.x[d] == "A")
  b2 <- which(combined.CLOZUK.PGC$A2.x[d] == "C")
  b3 <- which(combined.CLOZUK.PGC$A2.x[d] == "G")
  b4 <- which(combined.CLOZUK.PGC$A2.x[d] == "T")
  
  combined.CLOZUK.PGC$A2.x[d[b1]] <-"T"
  combined.CLOZUK.PGC$A2.x[d[b2]] <-"G"
  combined.CLOZUK.PGC$A2.x[d[b3]] <-"C"
  combined.CLOZUK.PGC$A2.x[d[b4]] <-"A"
  
  a <- which(combined.CLOZUK.PGC$A1.y==combined.CLOZUK.PGC$A2.x & combined.CLOZUK.PGC$A2.y==combined.CLOZUK.PGC$A1.x)
  # combined.CLOZUK.PGC$BP[a] <- (-combined.CLOZUK.PGC$BP[a])
  combined.CLOZUK.PGC$A1.x[a] <-as.character(combined.CLOZUK.PGC$A1.y[a])
  combined.CLOZUK.PGC$A2.x[a] <-as.character(combined.CLOZUK.PGC$A2.y[a])
  
  a <- which(combined.CLOZUK.PGC$A1.y != combined.CLOZUK.PGC$A1.x | combined.CLOZUK.PGC$A2.y!=combined.CLOZUK.PGC$A2.x)
  if (length(a)>0) combined.CLOZUK.PGC <- combined.CLOZUK.PGC[-a,]
  
}
#####

