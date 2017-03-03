# Checking the SNPs used in pathway analysis #

setwd("~/Documents/testing_PRS_chromosome_22/test_chr5/output/")

p.value.thresholds <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
x <- matrix()
abnormal_learning_memory_conditioning_score_files <- list(x,x,x,x,x,x,x,x,x,x,x,x)
row_number <- rep(0,12)
row_number2 <- rep(0,12) 


for (i in 1:length(p.value.thresholds)) {
current_matrix <- fread(paste0("./abnormal_learning_memory_conditioning/score/abnormal_learning_memory_conditioning_with_",p.value.thresholds[i],".score"))
assign(paste0("abnormal_learning_memory_conditioning_score_files",p.value.thresholds[i]), current_matrix, envir = .GlobalEnv)
row_number[i] <- nrow(current_matrix)
}

a <- merge(abnormal_learning_memory_conditioning_score_files0.01, abnormal_learning_memory_conditioning_score_files0.05, all = T, by = 'V1')
Shared_snps <- merge(abnormal_learning_memory_conditioning_score_files0.01, abnormal_learning_memory_conditioning_score_files0.05, all.x = T, all.y = F, by = 'V1')
unshared_SNPs <-a[is.na(V2.x),] 

for (i in 1:length(p.value.thresholds)) {
  current_matrix <- fread(paste0("./abnormal_long_term_potentiation/score/abnormal_long_term_potentiation_with_",p.value.thresholds[i],".score"))
  assign(paste0("abnormal_long_term_potentiation_score_files",p.value.thresholds[i]), current_matrix, envir = .GlobalEnv)
  row_number2[i] <- nrow(current_matrix)
}

a2 <- merge(abnormal_long_term_potentiation_score_files0.01, abnormal_long_term_potentiation_score_files0.05, all = T, by = 'V1')
Shared_snps2 <- merge(abnormal_long_term_potentiation_score_files0.01, abnormal_long_term_potentiation_score_files0.05, all.x = T, all.y = F, by = 'V1')
unshared_SNPs2 <-a2[is.na(V2.x),] 

sd(Shared_snps2$V3.y)
sd(unshared_SNPs2$V3.y)
sd(unshared_SNPs$V3.y)
sd(Shared_snps$V3.y)




