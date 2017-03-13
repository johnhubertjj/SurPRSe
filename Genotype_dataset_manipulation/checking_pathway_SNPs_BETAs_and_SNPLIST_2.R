# Checking the SNPs used in pathway analysis #

library(data.table)
setwd("~/Documents/testing_PRS_chromosome_22/test_chr5/output/")
Useful_pathways <- c("FMRP_targets", "abnormal_behavior", "abnormal_nervous_system_electrophysiology", "abnormal_learning_memory_conditioning", "abnormal_CNS_synaptic_transmission", "Cav2_channels", "abnormal_synaptic_transmission", "5HT_2C", "abnormal_long_term_potentiation", "abnormal_motor_capabilities_coordination_movement", "abnormal_behavioral_response_to_xenobiotic", "abnormal_associative_learning", "Lek2015_LoFintolerant_90", "BGS_top2_mean", "BGS_top2_max")
p.value.thresholds <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)


for (l in 1:length(Useful_pathways)) {  
x <- matrix()

abnormal_learning_memory_conditioning_score_files <- list(x,x,x,x,x,x,x,x,x,x,x,x)
row_number <- rep(0,12)
row_number2 <- rep(0,12) 
path_to_input_file <- paste0("./",Useful_pathways[l],"/score/")


for (i in 1:length(p.value.thresholds)) {
current_matrix <- fread(paste0(path_to_input_file, Useful_pathways[l],"_with_",p.value.thresholds[i],"_removing_SE_equaltoandmorethan_five.score"))
assign(paste0(Useful_pathways[l],"_score_files",p.value.thresholds[i]), current_matrix, envir = .GlobalEnv)
row_number[i] <- nrow(current_matrix)
}
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

a0.001 <- merge(`5HT_2C_score_files1e-04`, `5HT_2C_score_files0.001`, all = T, by = 'V1')
shared_snps0.001 <- merge(`5HT_2C_score_files1e-04`, `5HT_2C_score_files0.001`, all.x = T, all.y = F, by = 'V1')
unshared_snps0.001 <-a0.001[is.na(V2.x),] 

a0.01<- merge(`5HT_2C_score_files0.001`, `5HT_2C_score_files0.01`, all = T, by = 'V1')
shared_snps0.01 <- merge(`5HT_2C_score_files0.001`, `5HT_2C_score_files0.01`, all.x = T, all.y = F, by = 'V1')
unshared_snps0.01 <-a0.01[is.na(V2.x),] 

a0.05<- merge(`5HT_2C_score_files0.01`, `5HT_2C_score_files0.05`, all = T, by = 'V1')
shared_snps0.05 <- merge(`5HT_2C_score_files0.01`, `5HT_2C_score_files0.05`, all.x = T, all.y = F, by = 'V1')
unshared_snps0.05 <-a0.05[is.na(V2.x),] 

a0.1<- merge(`5HT_2C_score_files0.05`, `5HT_2C_score_files0.1`, all = T, by = 'V1')
shared_snps0.1 <- merge(`5HT_2C_score_files0.05`, `5HT_2C_score_files0.1`, all.x = T, all.y = F, by = 'V1')
unshared_snps0.1 <-a0.1[is.na(V2.x),] 

sd(shared_snps0.001$V3.y)
sd(unshared_snps0.001$V3.y)

mean(abs(shared_snps0.001$V3.y))
mean(abs(unshared_snps0.001$V3.y))

sd(shared_snps0.01$V3.y)
sd(unshared_snps0.01$V3.y)

mean(abs(shared_snps0.01$V3.y))
mean(abs(unshared_snps0.01$V3.y))

sd(shared_snps0.05$V3.y)
sd(unshared_snps0.05$V3.y)

mean(abs(shared_snps0.05$V3.y))
mean(abs(unshared_snps0.05$V3.y))

sd(shared_snps0.1$V3.y)
sd(unshared_snps0.1$V3.y)

mean(abs(shared_snps0.1$V3.y))
mean(abs(unshared_snps0.1$V3.y))



for (l in 1:length(Useful_pathways)) {  
  x <- matrix()
  
  abnormal_learning_memory_conditioning_score_files <- list(x,x,x,x,x,x,x,x,x,x,x,x)
  row_number <- rep(0,12)
  row_number2 <- rep(0,12) 
  path_to_input_file <- paste0("./",Useful_pathways[l],"/score/")
  
  
  for (i in 1:length(p.value.thresholds)) {
    current_matrix <- fread(paste0(path_to_input_file, Useful_pathways[l],"_with_",p.value.thresholds[i],".score"))
    assign(paste0(Useful_pathways[l],"_score_files",p.value.thresholds[i]), current_matrix, envir = .GlobalEnv)
    row_number[i] <- nrow(current_matrix)
  }
}