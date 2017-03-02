library(data.table)

abc <- fread("abnormal_associative_learning/score/abnormal_associative_learning_with_0.001_removing_SE_equaltoandmorethan_five_SNPS_in_question.txt")
abcd<- fread("abnormal_learning_memory_conditioning/score/abnormal_learning_memory_conditioning_with_0.001_removing_SE_equaltoandmorethan_five_SNPS_in_question.txt")

total2 <- abc[V11 >= 0.01 & V12 >= 0.01 & V13 >= 0.9]
total3 <- abcd[V11 > 0.01 & V12 > 0.01 & V13 > 0.9]

write.table(total3,file = "/Users/johnhubert/Dropbox/SNPs_which_have_high_SE_but_pass_MAF_INFO_abnormal_learning_memory_conditioning.txt", quote = F, row.names = F)
