#####################################
###########PRS PER PATHWAY###########
#####################################

setwd("~/Documents/testing_PRS_chromosome_22/test_chr5/output/")

library("data.table")
#SCORING#
Useful_pathways <- c("FMRP_targets", "abnormal_behavior", "abnormal_nervous_system_electrophysiology", "abnormal_learning_memory_conditioning", "abnormal_CNS_synaptic_transmission", "Cav2_channels", "abnormal_synaptic_transmission", "5HT_2C", "abnormal_long_term_potentiation", "abnormal_motor_capabilities_coordination_movement", "abnormal_behavioral_response_to_xenobiotic", "abnormal_associative_learning", "Lek2015_LoFintolerant_90", "BGS_top2_mean", "BGS_top2_max")
PGC_final <- fread("combined_PGC_table_with_CHR.POS_identifiers.txt")
p.value.thresholds <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

for (i in 1:length(Useful_pathways)) {
  Maindir <- paste0("~/Documents/testing_PRS_chromosome_22/test_chr5/output/", Useful_pathways[i])
  scoredir <- "score"
  profiledir <- "Profile"
  if (file.exists(scoredir) == FALSE){
    dir.create(file.path(Maindir, scoredir))
  }
  if (file.exists(profiledir) == FALSE){
    dir.create(file.path(Maindir, profiledir))
  }
  
  scoring_output_file <- paste0(Useful_pathways[i],"/scoring_PGC_CLOZUK_pathway", Useful_pathways[i])
  Current_pathway <- fread(paste0(Useful_pathways[i],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[i],".bim")) 
  path_to_pathway_plink_file <- paste0(Useful_pathways[i],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[i])                         
  setnames(Current_pathway,c("CHR","SNP","GD","BP","A1","A2"))                        
  combined.CLOZUK.PGC.clumped.Genomic.SNPs <- merge(Current_pathway, PGC_final, by.x="SNP", by.y="SNP", all=F, sort=F)
  combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1.y <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A1.y)
  combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2.y <- toupper(combined.CLOZUK.PGC.clumped.Genomic.SNPs$A2.y)
  
  ## Change so that the SNPs with SE >= 5 are removed
  SNPs_removed <- combined.CLOZUK.PGC.clumped.Genomic.SNPs [SE >= 5]
  combined.CLOZUK.PGC.clumped.Genomic.SNPs <- combined.CLOZUK.PGC.clumped.Genomic.SNPs[!SE >= 5]
  
  ## check that it is merging properly here (after analysis is run)
  
  for (w in 1:length(p.value.thresholds)) {
    
    a <- copy(combined.CLOZUK.PGC.clumped.Genomic.SNPs)    
    SNPs <- a[, .I[which(P <= p.value.thresholds[w])]]
    
    if (length(SNPs) != 0){
      a <- a[SNPs, .(SNP,A1.y,BETA)]
      
      filename <- paste0(Useful_pathways[i],'/score/', Useful_pathways[i],'_with_', p.value.thresholds[w],"_removing_SE_equaltoandmorethan_five.score")
      filename_SNPs_removed <- paste0(Useful_pathways[i],'/score/', Useful_pathways[i],'_with_', p.value.thresholds[w],"_removing_SE_equaltoandmorethan_five_SNPS_in_question.txt")
      
      write.table(file = filename, a, row.names = F, col.names = F, quote = F, sep="\t")
      write.table(file = filename_SNPs_removed, SNPs_removed, row.names = F, col.names = F, quote = F, sep="\t")
      
      #      if (Useful_pathways[i] == "abnormal_learning|memory|conditioning" | Useful_pathways[i] == "abnormal_motor_capabilities|coordination|movement") {
      #        Useful_pathways[i] <- "abnormal_learning\\|memory\\|conditioning"
      #        filename <- paste0(Useful_pathways[i],'/score/', Useful_pathways[i],'_with_', p.value.thresholds[w],".score")
      #        path_to_pathway_plink_file <- paste0(Useful_pathways[i],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[i])                         
      #        }
      
      #      if (Useful_pathways[i] == "abnormal_motor_capabilities|coordination|movement") {
      #        Useful_pathways[i] <-"abnormal_learning_motor_capabilities\|coordination\|movement"
      #        filename <- paste0(Useful_pathways[i],'/score/', Useful_pathways[i],'_with_', p.value.thresholds[w],".score")
      #        path_to_pathway_plink_file <- paste0(Useful_pathways[i],"/pathways_CLOZUK_GWAS_BGE_CLUMPED_",Useful_pathways[i])                         
      #        
      #      }
      
      command <- paste0('/Users/johnhubert/Documents/plink_mac//plink --bfile ',path_to_pathway_plink_file, ' --score ', filename, " --out ", path_to_pathway_plink_file, "_", p.value.thresholds[w], "_removing_SE_equaltoandmorethan_five")
      system(command)
      rm(a)
    }else{
      next()
    }
  }
}

write.table(Genes_used, file = "Index_of_genes_and_pval_1.txt", quote = F, col.names = T, row.names = F)
#score

for (i in 1:(length(p.value.thresholds)))
{
  # format the p.values so that exponential notation is NOT used, could result in confusion in BASH scripting
  scoring.output.filename <- paste0(scoring_output_file, "_", format(p.value.thresholds[i], scientific =  F), ".score")
  SNPs.lower.than.threshold <- which(combined.CLOZUK.PGC.clumped.Genomic.SNPs$P.x <= p.value.thresholds[i]) 
  
  if (length(SNPs.lower.than.threshold ) > 0) {
    tmp <- combined.CLOZUK.PGC.clumped.Genomic.SNPs[SNPs.lower.than.threshold, ]
    final.output <- tmp[, c("SNP","A1","BETA"), with = F]
  }
  
  if (nrow(out)>0) {
    write.table(file = scoring.output.filename, final.output, row.names = F, col.names = F, quote = F, sep="\t")
  }
}
## END script!
#######################################