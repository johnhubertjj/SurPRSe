# No longer required, but may be useful to print off plots for results 
```{r}
# for (t in 1:4){
Final_analysis <- list(
  current_table_extended <- Full_data[grepl("^extended",Full_data$score)],
  current_table_normal <- Full_data[grepl("^normal",Full_data$score)],
  current_table_extended_All.genome <- Full_data[grepl("^extended",Full_data$score) | grepl("All.genome", Full_data$score)],
  current_table_normal_All.genome <- Full_data[grepl("^normal", Full_data$score) | grepl("All.genome", Full_data$score)]
)

#current_table_extended <- current_table_extended[order(current_table$p)]
#current_table_normal <- current_table_normal[order(current_table_normal$p)]
#current_table_extended_All.genome <- current_table_extended_All.genome[order(current_table_extended_All.genome$p)]
#current_table_normal_All.genome <- current_table_normal_All.genome[order(current_table_normal_All.genome$p)]



#Full_analysis_list <- list( analysis1, analysis2, analysis3, analysis4)
#cut_down_analysis_list <- Full_analysis_list

Analysis_vector <- c("normal vs whole genome genic", "extended vs whole genome genic", "normal vs whole genome all", "extended vs whole genome all")

samples2 <- c("Everything", "broad", "narrow", "All_minus_narrow")

for (t in 3:4){
  
  current_analysis <- Final_analysis[[t]]
  
  if (t > 2){
    to_remove <- current_analysis[,.I[grepl("genic.genome",current_analysis$score)]]
    current_analysis <- current_analysis[-to_remove]
  } 
  setkey(current_analysis, "samples.i.")
  cutdown_tables_current <- current_analysis[order(p)]
  cutdown_tables_current <- cutdown_tables_current[,.SD[1:10], by = samples.i.]
  
  cutdown_tables_current %>% 
    mutate_if(is.factor, as.character) -> cutdown_tables_current
  
  print(
    kable(cutdown_tables_current[,c(1,2,4:ncol(cutdown_tables_current))], format = "html",caption = paste0("Cognition;",Analysis_vector[t],";across different DSM ranges")) %>%
      kable_styling("striped", full_width = F, position = "left",bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
      group_rows(cutdown_tables_current[1,"samples.i."], 1, 10) %>%
      group_rows(cutdown_tables_current[11,"samples.i."], 11, 20) %>%
      group_rows(cutdown_tables_current[21,"samples.i."], 21, 30) %>%
      group_rows(cutdown_tables_current[31,"samples.i."], 31, 40 ) %>%  
      add_footnote(c("Some footnote", "Some other footnote"))
  )
}

for (t in 1:3){
  kable(all_TRS_PRS_results[[t]][1:20,], caption = "TRS tables across different DSM ranges")
}
```

### Plotting
```{r echo = F,message = F, results=F, warning=F}

Final_analysis_gglist_pval <- list()
Final_analysis_beta <- list()
Final_analysis_r2_pathways <- list()

gglist_pathways <- list()
gglist_beta_pathways <- list()
gglist_r2_dir_pathways <- list()

for (t in 1:4){
  
  current_table <- Final_analysis[[t]]
  
  # Pathways
  str(current_table)
  
  COGS_only_0.5 <- current_table[grepl("_0.5",current_table$score)] 
  COGS_only_0.05 <- current_table[grepl("_0.05",current_table$score)] 
  
  #COGS_only_0.05 <- COGS_only_0.05[!grepl("abnormal_grooming_behavior", COGS_only_0.05$score)]
  #COGS_only_0.5 <- COGS_only_0.5[!grepl("abnormal_grooming_behavior", COGS_only_0.5$score)]
  
  #COGS_only_0.05 <- COGS_only_0.05[!grepl("Lek2015_LoFintolerant_90", COGS_only_0.05$score)]
  #COGS_only_0.5 <- COGS_only_0.5[!grepl("Lek2015_LoFintolerant_90", COGS_only_0.5$score)]
  
  split_data_frame_0.05 <- split(COGS_only_0.05, by = "test")
  split_data_frame_0.5 <- split(COGS_only_0.5, by = "test")
  
  #setkey(split_data_frame_0.05,p)
  #setkey(split_data_frame_0.5,p)
  
  #Whole brain
  #ALSPAC_all_brain_0.05 <- ALSPAC_full_brain_regions[grepl("_0.05",ALSPAC_full_brain_regions$score)]
  #ALSPAC_all_brain_0.5 <- ALSPAC_full_brain_regions[grepl("_0.5",ALSPAC_full_brain_regions$score)] 
  
  
  
  for (i in 1:length(split_data_frame_0.05)){
    
    
    # Obtaining the log p-values
    top_five_pathways <- split_data_frame_0.05[[i]][order(split_data_frame_0.05[[i]]$score)]
    top_five_pathways$P_threshold <- 0.05
    
    # Alter to find the whole genome row from the pathway rows
    top_five_pathways_0.05 <- top_five_pathways[,.I[grepl("genic.genome",top_five_pathways$score) | grepl("All.genome", top_five_pathways$score)]]
    Type_vector <- rep("Pathway",nrow(top_five_pathways))
    Type_vector[top_five_pathways_0.05] <- "Whole_genome"
    top_five_pathways$Type <- Type_vector
    
    # By Factoring the p-values, you can get them to trend decreasingly without affecting the whole genome p-value
    top_five_pathways$p <- as.numeric(top_five_pathways$p)
    top_five_pathways$score <- factor(top_five_pathways$score, levels = top_five_pathways$score[order(top_five_pathways$score, top_five_pathways$Type)])
    
    
    # Do the same for the 0.5 p-thresholds
    top_five_pathways2 <- split_data_frame_0.5[[i]][order(split_data_frame_0.5[[i]]$score)]
    top_five_pathways2$P_threshold <- 0.5
    
    top_five_pathways_0.5 <- top_five_pathways2[,.I[grepl("genic.genome",top_five_pathways2$score) | grepl("All.genome", top_five_pathways2$score)]]
    Type_vector <- rep("Pathway",nrow(top_five_pathways2))
    Type_vector[top_five_pathways_0.5] <- "Whole_genome"
    top_five_pathways2$Type <- Type_vector
    
    
    # By Factoring the p-values, you can get them to trend decreasingly without affecting the whole genome p-value
    top_five_pathways2$p <- as.numeric(top_five_pathways2$p)
    top_five_pathways2$score <- factor(top_five_pathways2$score, levels = top_five_pathways2$score[order(top_five_pathways2$score, top_five_pathways2$Type)])
    
    
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
    whole_genome_positions <- grep(x = test1, pattern = "SCORE_whole_genome_\\d.\\d+",perl = T)
    
    if (t > 2){
      whole_genome_genic_positions <- grep(x = test1, pattern = "genic.genome_SCORE_whole_genome_\\d.\\d+",perl = T)
      whole_genome_all_genome_positions <- grep(x = test1, pattern = "All.genome_SCORE_whole_genome_\\d.\\d+",perl = T)
    }
    
    alterations <- test1
    alterations[-whole_genome_positions] <- str_replace(string = test1[-whole_genome_positions], pattern = "SCORE_(.*)_\\d.\\d+",replacement = "\\1")
    alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "long_term_potentiation",replacement = "LTP")
    alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "action_potential",replacement = "AP")
    alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "depolarization",replacement = "DP")
    alterations[whole_genome_positions] <- rep(x = "Whole Genome PRS", 2)
    
    if (t > 2){
      alterations[whole_genome_genic_positions] <- rep(x = "Whole Genome PRS GENIC", 2)
      alterations[whole_genome_all_genome_positions] <- rep(x = "Whole Genome PRS ALL", 2)
    }
    
    test2 <- vector(length = length(test1))
    names(alterations) <- test1
    
    p <- p + scale_x_discrete(labels= alterations)
    #p <- p + scale_y
    p <- p + facet_grid(. ~ P_threshold,scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    p <- p + geom_hline(aes(yintercept=1.30103), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + scale_fill_brewer(palette = "Paired")
    p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggtitle(top_five_pathways_all$.id[1])
    #p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5,vjust = 0.5,angle = 90))
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = expression(-'log'[10]*'(p)'))
    p <- p + xlab(label = "Polygenic risk score")
    
    gglist_pathways[[i]] <- p
    
    p <- ggplot(top_five_pathways_all, aes(x=score, y=estimate, fill = Type, group=P_threshold))
    
    p <- p +
      geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25) +
      geom_point(aes(colour = Type))
    
    test1 <- levels(p$data$score)
    whole_genome_positions <- grep(x = test1,pattern = "SCORE_whole_genome_\\d.\\d+",perl = T)
    
    if (t > 2){
      whole_genome_genic_positions <- grep(x = test1, pattern = "genic.genome_SCORE_whole_genome_\\d.\\d+",perl = T)
      whole_genome_all_genome_positions <- grep(x = test1, pattern = "All.genome_SCORE_whole_genome_\\d.\\d+",perl = T)
    }
    
    
    alterations <- test1
    alterations[-whole_genome_positions] <- str_replace(string = test1[-whole_genome_positions], pattern = "SCORE_(.*)_\\d.\\d+",replacement = "\\1")
    alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "long_term_potentiation",replacement = "LTP")
    alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "action_potential",replacement = "AP")
    alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "depolarization",replacement = "DP")
    alterations[whole_genome_positions] <- rep(x = "Whole Genome PRS", 2)
    
    if (t > 2){
      alterations[whole_genome_genic_positions] <- rep(x = "Whole Genome PRS GENIC", 2)
      alterations[whole_genome_all_genome_positions] <- rep(x = "Whole Genome PRS ALL", 2)
    }
    
    test2 <- vector(length = length(test1))
    names(alterations) <- test1
    
    p <- p + scale_x_discrete(labels= alterations)
    p <- p + facet_grid(. ~ P_threshold, scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + scale_fill_brewer(palette = "Paired")
    p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggtitle(top_five_pathways_all$.id[1])
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = "BETA")
    p <- p + xlab(label = "Polygenic risk score")
    
    gglist_beta_pathways[[i]] <- p
    
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
    whole_genome_positions <- grep(x = test1,pattern = "SCORE_whole_genome_\\d.\\d+",perl = T)
    alterations <- test1
    
    if (t > 2){
      whole_genome_genic_positions <- grep(x = test1, pattern = "genic.genome_SCORE_whole_genome_\\d.\\d+",perl = T)
      whole_genome_all_genome_positions <- grep(x = test1, pattern = "All.genome_SCORE_whole_genome_\\d.\\d+",perl = T)
    }
    
    
    alterations[-whole_genome_positions] <- str_replace(string = test1[-whole_genome_positions], pattern = "SCORE_(.*)_\\d.\\d+",replacement = "\\1")
    alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "long_term_potentiation",replacement = "LTP")
    alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "action_potential",replacement = "AP")
    alterations[-whole_genome_positions] <- str_replace(string = alterations[-whole_genome_positions], pattern = "depolarization",replacement = "DP")
    alterations[whole_genome_positions] <- rep(x = "Whole Genome PRS", 2)
    test2 <- vector(length = length(test1))
    names(alterations) <- test1
    
    if (t > 2){
      alterations[whole_genome_genic_positions] <- rep(x = "Whole Genome PRS GENIC", 2)
      alterations[whole_genome_all_genome_positions] <- rep(x = "Whole Genome PRS ALL", 2)
    }
    
    p <- p + scale_x_discrete(labels= alterations)
    
    p <- p + scale_y_continuous(expand = expand_scale(mult = c(0.2,.6)))
    p <- p + facet_grid(. ~ P_threshold, scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    p <- p + theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1,vjust = 0.5))
    p <- p + ggtitle(top_five_pathways_all$.id[1])
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = "R2_dir (%)")
    p <- p + xlab(label = "Polygenic risk score")
    
    gglist_r2_dir_pathways[[i]] <- p
    
    
  }
  
  Final_analysis_gglist_pval[[t]] <- gglist_pathways
  Final_analysis_beta[[t]] <- gglist_beta_pathways
  Final_analysis_r2_pathways[[t]] <- gglist_r2_dir_pathways
  
}





```
#### Save the plots to File

```{r, message=F, include=F}
pdf("~/Documents/COGS_pathways_results/COGS_CLOZUK_pathway_whole_genome_genic_comparison_R2_dir_12_pathways_EXTENDED_REGIONS.pdf", onefile = TRUE)
invisible(lapply(Final_analysis_r2_pathways[[3]],print))
dev.off()

pdf("~/Documents/COGS_pathways_results/COGS_CLOZUK_pathway_whole_genome_genic_comparison_R2_dir_12_pathways_NORMAL_REGIONS.pdf", onefile = TRUE)
invisible(lapply(Final_analysis_r2_pathways[[4]],print))
dev.off()


pdf("~/Documents/COGS_pathways_results/COGS_CLOZUK_pathway_whole_genome_genic_comparison_beta_12_pathways_EXTENDED_REGIONS.pdf", onefile = TRUE)
invisible(lapply(Final_analysis_beta[[3]],print))
dev.off()

pdf("~/Documents/COGS_pathways_results/COGS_CLOZUK_pathway_whole_genome_genic_comparison_beta_12_pathways_NORMAL_REGIONS.pdf", onefile = TRUE)
invisible(lapply(Final_analysis_beta[[4]],print))
dev.off()

# NOT RUN {
x <- knitr::kable(head(mtcars), "html")
# Add a row of header with 3 columns on the top of the table. The column
# span for the 2nd and 3rd one are 5 & 6.
add_header_above(x, c(" ", "Group 1" = 5, "Group 2" = 6))

```
