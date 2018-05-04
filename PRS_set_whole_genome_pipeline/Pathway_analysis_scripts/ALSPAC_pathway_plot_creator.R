load("~/Documents/ALSPAC_paper_workspaces/Brain_regions_common_risk.RData")
library(data.table)
# required packages -- will only install if they are not already installed
list.of.packages <- c("plyr", "stringr", "dplyr", "tidyr", "reshape2", "ggplot2", "scales", "data.table", "plotly", "devtools", "gridExtra", "cowplot", "repr", "knitr", "kableExtra", "IRdisplay")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, "/home/c1020109/jupyter-notebooks/lib/R/library")


#library(ggplot2) 

# loads the required packages
lapply(list.of.packages, require, character.only = TRUE)

ALSPAC_all_pathways <- fread("~/Documents/ALSPAC_paper_workspaces/CLOZUK_PGC2noclo_Biobank_Selected_SCZ_pathways.txt")


ALSPAC_full_brain_regions <- fread("~/Documents/ALSPAC_paper_workspaces/CLOZUK_PGC2noclo_biobank_whole_genome.txt")
testa <- ALSPAC_full_brain_regions


# Pathways
str(ALSPAC_all_pathways)

ALSPAC_only_0.5 <- ALSPAC_all_pathways[grepl("_0.5",ALSPAC_all_pathways$score)] 
ALSPAC_only_0.05 <- ALSPAC_all_pathways[grepl("_0.05",ALSPAC_all_pathways$score)] 

ALSPAC_only_0.05 <- ALSPAC_only_0.05[!grepl("abnormal_grooming_behavior", ALSPAC_only_0.05$score)]
ALSPAC_only_0.5 <- ALSPAC_only_0.5[!grepl("abnormal_grooming_behavior", ALSPAC_only_0.5$score)]

ALSPAC_only_0.05 <- ALSPAC_only_0.05[!grepl("Lek2015_LoFintolerant_90", ALSPAC_only_0.05$score)]
ALSPAC_only_0.5 <- ALSPAC_only_0.5[!grepl("Lek2015_LoFintolerant_90", ALSPAC_only_0.5$score)]

ALSPAC_only_0.05 <- ALSPAC_only_0.05[!grepl("grey", ALSPAC_only_0.05$test)]
ALSPAC_only_0.5 <- ALSPAC_only_0.5[!grepl("grey", ALSPAC_only_0.5$test)]

split_data_frame_0.05 <- split(ALSPAC_only_0.05, by = "test")
split_data_frame_0.5 <- split(ALSPAC_only_0.5, by = "test")
#setkey(split_data_frame_0.05,p)
#setkey(split_data_frame_0.5,p)

#Whole brain
ALSPAC_all_brain_0.05 <- ALSPAC_full_brain_regions[grepl("_0.05",ALSPAC_full_brain_regions$score)]
ALSPAC_all_brain_0.5 <- ALSPAC_full_brain_regions[grepl("_0.5",ALSPAC_full_brain_regions$score)] 
gglist_pathways_biobank <- list()
gglist_beta_pathways_biobank <- list()
gglist_r2_dir_pathways_biobank <- list()

for (i in 1:length(split_data_frame_0.05)){
  
  # Obtaining the log p-values
  top_five_pathways <- split_data_frame_0.05[[i]][order(split_data_frame_0.05[[i]]$score)]
  #top_five_pathways <- top_five_pathways[1:5]
  whole_brain_p_value_0.05 <- ALSPAC_all_brain_0.05[which(top_five_pathways$test[1] == ALSPAC_all_brain_0.05)]
  top_five_pathways <- rbind(top_five_pathways,whole_brain_p_value_0.05)
  top_five_pathways$P_threshold <- 0.05
  top_five_pathways$Type <- c(rep("Pathway",9),"Whole_genome")
  top_five_pathways_0.05 <- which(top_five_pathways$score == "SCORE_0.05")
  
  # By Factoring the p-values, you can get them to trend decreasingly without affecting the whole genome p-value
  top_five_pathways$p <- as.numeric(top_five_pathways$p)
  top_five_pathways$score <- factor(top_five_pathways$score, levels = top_five_pathways$score[c(top_five_pathways_0.05, order(top_five_pathways$score[-top_five_pathways_0.05]))])
  
  
  # Do the same for the 0.5 p-thresholds
  top_five_pathways2 <- split_data_frame_0.5[[i]][order(split_data_frame_0.5[[i]]$score)]
  #top_five_pathways2 <- top_five_pathways2[1:5]
  whole_brain_p_value_0.5 <- ALSPAC_all_brain_0.5[which(top_five_pathways2$test[1] == ALSPAC_all_brain_0.5)]
  top_five_pathways2 <- rbind(top_five_pathways2,whole_brain_p_value_0.5)
  top_five_pathways2$P_threshold <- 0.5
  top_five_pathways2$Type <- c(rep("Pathway",9),"Whole_genome")
  
  top_five_pathways_0.5 <- which(top_five_pathways2$score == "SCORE_0.5")
  top_five_pathways2$p <- as.numeric(top_five_pathways2$p)
  top_five_pathways2$score <- factor(top_five_pathways2$score, levels = top_five_pathways2$score[c(top_five_pathways_0.5, order(top_five_pathways2$score[-top_five_pathways_0.5]))])
  
  
  top_five_pathways_all <- rbind(top_five_pathways,top_five_pathways2)
  top_five_pathways_all$p <- as.numeric(top_five_pathways_all$p)
  top_five_pathways_all$logp <- -log10(top_five_pathways_all$p)
  
  top_five_pathways_all$estimate <- as.numeric(top_five_pathways_all$estimate)
  top_five_pathways_all$upper <- as.numeric(top_five_pathways_all$upper)
  top_five_pathways_all$lower <- as.numeric(top_five_pathways_all$lower)
  top_five_pathways_all$SE <- as.numeric(top_five_pathways_all$SE)
  
  top_five_pathways_all$SE_higher <- top_five_pathways_all$estimate + top_five_pathways_all$SE
  top_five_pathways_all$SE_lower <- top_five_pathways_all$estimate - top_five_pathways_all$SE
  
  top_five_pathways_all$r2_dir <- 100 * (as.numeric(top_five_pathways_all$r.squared) *
                                           (sign(as.numeric(top_five_pathways_all$estimate))))
  
  top_five_pathways_all$p_value_text <- paste("p =", scientific(top_five_pathways_all$p, digits = 2), sep = " ")
  assign(paste0(top_five_pathways_all$test[1],"_top_five_pathways_all_thresholds"), top_five_pathways_all, envir = .GlobalEnv)
  
  p <- ggplot(top_five_pathways_all, aes(x=score, y=logp, fill = Type, group=P_threshold))
  
  p <- p +
    geom_point(aes(colour = Type))
  
  test1 <- levels(p$data$score)
  whole_genome_positions <- grep(x = test1,pattern = "SCORE_\\d.\\d+",perl = T)
  alterations <- test1
  alterations[-whole_genome_positions] <- str_replace(string = test1[-whole_genome_positions], pattern = "SCORE_(.*)_\\d.\\d+",replacement = "\\1")
  alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "long_term_potentiation",replacement = "LTP")
  alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "action_potential",replacement = "AP")
  alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "depolarization",replacement = "DP")
  alterations[whole_genome_positions] <- rep(x = "Whole Genome PRS", 2)
  test2 <- vector(length = length(test1))
  names(alterations) <- test1
  
  p <- p + scale_x_discrete(labels= alterations)
  #p <- p + scale_y
  p <- p + facet_grid(. ~ P_threshold,scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  p <- p + geom_hline(aes(yintercept=1.30103), colour = "red", linetype= "solid", alpha = 0.25)
  p <- p + scale_fill_brewer(palette = "Paired")
  p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
  p <- p + ggtitle(top_five_pathways_all$test[1])
  #p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5,vjust = 0.5,angle = 90))
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  p <- p + ylab(label = expression(-'log'[10]*'(p)'))
  p <- p + xlab(label = "Polygenic risk score")

  gglist_pathways_biobank[[i]] <- p
  
  p <- ggplot(top_five_pathways_all, aes(x=score, y=estimate, fill = Type, group=P_threshold))
  
  p <- p +
    geom_errorbar(aes(ymin = SE_higher, ymax = SE_lower), position = "dodge", width = 0.25) +
    geom_point(aes(colour = Type))
  
  test1 <- levels(p$data$score)
  whole_genome_positions <- grep(x = test1,pattern = "SCORE_\\d.\\d+",perl = T)
  alterations <- test1
  alterations[-whole_genome_positions] <- str_replace(string = test1[-whole_genome_positions], pattern = "SCORE_(.*)_\\d.\\d+",replacement = "\\1")
  alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "long_term_potentiation",replacement = "LTP")
  alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "action_potential",replacement = "AP")
  alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "depolarization",replacement = "DP")
  alterations[whole_genome_positions] <- rep(x = "Whole Genome PRS", 2)
  test2 <- vector(length = length(test1))
  names(alterations) <- test1
  
  p <- p + scale_x_discrete(labels= alterations)
  p <- p + facet_grid(. ~ P_threshold, scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
  p <- p + scale_fill_brewer(palette = "Paired")
  p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
  p <- p + ggtitle(top_five_pathways_all$test[1])
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  p <- p + ylab(label = "BETA")
  p <- p + xlab(label = "Polygenic risk score")
  
  gglist_beta_pathways_biobank[[i]] <- p
  
  p <- ggplot(top_five_pathways_all, aes(x=score, y=r2_dir, fill = Type, group=P_threshold))
  p <- p +
    geom_bar(stat = "identity", aes(colour = Type), position = "dodge") +
    geom_text(data=subset(top_five_pathways_all, p < 0.05),
              aes(x=score,y=r2_dir,label=p_value_text, hjust=ifelse(sign(r2_dir)>0, 0, 0)), angle = 90, position = position_dodge(width = 1), size = 2.9)

 #Problem with labels with a workaround
  # I use the score column in the format of factors and reference each relevant dataset for ggplot.
  # However this relies on having 0.05 and 0.5 in the value name.
  # scale_x_discrete accepts functions, but I also need to convert SCORE_0.05 and Score_0.5 into a "Whole_genome_PRS" which is almost impossible to write"
  # However as the labels function accepts key:value pairs, I wrote a vector in R that maps the original names of the pathways to "human readable" format using names function in R
  # This should work for most instances
  
  test1 <- levels(p$data$score)
  whole_genome_positions <- grep(x = test1,pattern = "SCORE_\\d.\\d+",perl = T)
  alterations <- test1
  alterations[-whole_genome_positions] <- str_replace(string = test1[-whole_genome_positions], pattern = "SCORE_(.*)_\\d.\\d+",replacement = "\\1")
  alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "long_term_potentiation",replacement = "LTP")
  alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "action_potential",replacement = "AP")
  alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "depolarization",replacement = "DP")
  alterations[whole_genome_positions] <- rep(x = "Whole Genome PRS", 2)
  test2 <- vector(length = length(test1))
  names(alterations) <- test1
  
  p <- p + scale_x_discrete(labels= alterations)

  p <- p + scale_y_continuous(expand = expand_scale(mult = c(0.2,.6)))
  p <- p + facet_grid(. ~ P_threshold, scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  p <- p + theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1,vjust = 0.5))
  p <- p + ggtitle(top_five_pathways_all$test[1])
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  p <- p + ylab(label = "R2_dir (%)")
  p <- p + xlab(label = "Polygenic risk score")
  
  gglist_r2_dir_pathways_biobank[[i]] <- p
  
}

pdf("/neurocluster/filesync/c1020109/Biobank_results/Pathway_results/Biobank_CLOZUK_pathway_whole_genome_genic_comparison_R2_dir_9_pathways.pdf", onefile = TRUE)

pdf("~/Documents/ALSPAC_paper_workspaces/Biobank_CLOZUK_pathway_whole_genome_genic_comparison_R2_dir_9_pathways.pdf", onefile = TRUE)
invisible(lapply(gglist_r2_dir_pathways_biobank,print))
dev.off()

