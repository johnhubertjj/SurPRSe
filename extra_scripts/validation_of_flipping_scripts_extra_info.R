test_validatation <- fread("validation.bim")
test_validatation$V5 <- toupper(test_validatation$V5)
test_validatation$V6 <- toupper(test_validatation$V6)

a <- which(test_validatation$V5 == "A" & test_validatation$V6 == "T")
if (length(a) != 0){
  test_validatation[-a]
}else{
  warning("there are no SNPs with A-T in dataset")
}
a <- which(test_validatation$V5 == "T" & test_validatation$V6 == "A")
if (length(a) != 0){
  test_validatation[-a]
}else{
  warning("there are no SNPs with T-A in dataset")
}

a <- which(test_validatation$V5 == "C" & test_validatation$V6 == "G")

if (length(a) != 0){
  test_validatation[-a]
}else{
  warning("there are no SNPs with C-G in dataset")
}

a <- which(test_validatation$V5 == "G" & test_validatation$V6 == "C")

if (length(a) != 0){
  test_validatation[-a]
}else{
  warning("there are no SNPs with G-C in dataset")
}

nrow(test_validatation)
colnames(test_validatation) <- c("CHR", "SNP", "GD", "BP", "A1", "A2")
data <- fread("training.assoc.logistic.snps")
merging_data <- merge(test_validatation,data, by = "SNP", all = F, sort = F)
nrow(merging_data)
write.table(file = "common.snps", merging_data$SNP, row.names = F, quote = F, col.names = F, sep = '\t')
