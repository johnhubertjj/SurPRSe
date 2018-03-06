library(data.table)

#ALSPAC_all_pathways <- fread(file = "CLOZUK_PGC2noclo_all_pocklington_inclu_LoF_specific_scores.txt")
ALSPAC_all_pathways <- fread("~/Dropbox/Untitled Folder/CLOZUK_PGC2noclo_Selected_SCZ_pathways.txt")

#ALSPAC_full_brain_regions <- results_CLOZUK_PGC2noclo_more_sig_thresh
ALSPAC_full_brain_regions <- results_CLOZUK_PGC2noclo
testa <- ALSPAC_full_brain_regions
ALSPAC_full_brain_regions <- as.data.table(ALSPAC_full_brain_regions)

# Pathways
str(ALSPAC_all_pathways)

ALSPAC_only_0.5 <- ALSPAC_all_pathways[grepl("_0.5",ALSPAC_all_pathways$score)] 
ALSPAC_only_0.05 <- ALSPAC_all_pathways[grepl("_0.05",ALSPAC_all_pathways$score)] 

ALSPAC_only_0.05 <- ALSPAC_only_0.05[!grepl("abnormal_grooming_behavior", ALSPAC_only_0.05$score)]
ALSPAC_only_0.5 <- ALSPAC_only_0.5[!grepl("abnormal_grooming_behavior", ALSPAC_only_0.5$score)]

ALSPAC_only_0.05 <- ALSPAC_only_0.05[!grepl("Lek2015_LoFintolerant_90", ALSPAC_only_0.05$score)]
ALSPAC_only_0.5 <- ALSPAC_only_0.5[!grepl("Lek2015_LoFintolerant_90", ALSPAC_only_0.5$score)]

split_data_frame_0.05 <- split(ALSPAC_only_0.05, by = "test")
split_data_frame_0.5 <- split(ALSPAC_only_0.5, by = "test")
setkey(split_data_frame_0.05,p)
setkey(split_data_frame_0.5,p)

#Whole brain
ALSPAC_all_brain_0.05 <- ALSPAC_full_brain_regions[grepl("_0.05",ALSPAC_full_brain_regions$score)]
ALSPAC_all_brain_0.5 <- ALSPAC_full_brain_regions[grepl("_0.5",ALSPAC_full_brain_regions$score)] 
gglist_pathways <- list()
gglist_beta_pathways <- list()
gglist_r2_dir_pathways <- list()

for (i in 1:length(split_data_frame_0.05)){
  
  # Obtaining the log p-values
  top_five_pathways <- split_data_frame_0.05[[i]][order(split_data_frame_0.05[[i]]$p)]
  top_five_pathways <- top_five_pathways[1:5]
  whole_brain_p_value_0.05 <- ALSPAC_all_brain_0.05[which(top_five_pathways$test[1] == ALSPAC_all_brain_0.05)]
  top_five_pathways <- rbind(top_five_pathways,whole_brain_p_value_0.05)
  top_five_pathways$P_threshold <- 0.05
  top_five_pathways$Type <- c(rep("Pathway",5),"Whole_genome")
  top_five_pathways_0.05 <- which(top_five_pathways$score == "SCORE_0.05")
  
  # By Factoring the p-values, you can get them to trend decreasingly without affecting the whole genome p-value
  top_five_pathways$p <- as.numeric(top_five_pathways$p)
  top_five_pathways$score <- factor(top_five_pathways$score, levels = top_five_pathways$score[c(top_five_pathways_0.05, order(top_five_pathways$p[-top_five_pathways_0.05]))])
  
  
  # Do the same for the 0.5 p-thresholds
  top_five_pathways2 <- split_data_frame_0.5[[i]][order(split_data_frame_0.5[[i]]$p)]
  top_five_pathways2 <- top_five_pathways2[1:5]
  whole_brain_p_value_0.5 <- ALSPAC_all_brain_0.5[which(top_five_pathways2$test[1] == ALSPAC_all_brain_0.5)]
  top_five_pathways2 <- rbind(top_five_pathways2,whole_brain_p_value_0.5)
  top_five_pathways2$P_threshold <- 0.5
  top_five_pathways2$Type <- c(rep("Pathway",5),"Whole_genome")
  
  top_five_pathways_0.5 <- which(top_five_pathways2$score == "SCORE_0.5")
  top_five_pathways2$p <- as.numeric(top_five_pathways2$p)
  top_five_pathways2$score <- factor(top_five_pathways2$score, levels = top_five_pathways2$score[c(top_five_pathways_0.5, order(top_five_pathways2$p[-top_five_pathways_0.5]))])
  
  
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
  
  p <- p + scale_x_discrete(labels=c("SCORE_0.05" = "Whole Genome PRS", "SCORE_0.5" = "Whole Genome PRS"))
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
  p
  gglist_pathways[[i]] <- p
  
  p <- ggplot(top_five_pathways_all, aes(x=score, y=estimate, fill = Type, group=P_threshold))
  
  p <- p +
    geom_errorbar(aes(ymin = SE_higher, ymax = SE_lower), position = "dodge", width = 0.25) +
    geom_point(aes(colour = Type))
  
  p <- p + scale_x_discrete(labels=c("SCORE_0.05" = "Whole Genome PRS", "SCORE_0.5" = "Whole Genome PRS"))
  p <- p + facet_grid(. ~ P_threshold, scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
  p <- p + scale_fill_brewer(palette = "Paired")
  p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
  p <- p + ggtitle(top_five_pathways_all$test[1])
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  p <- p + ylab(label = "BETA")
  p <- p + xlab(label = "Polygenic risk score")
  
  gglist_beta_pathways[[i]] <- p

  p <- ggplot(top_five_pathways_all, aes(x=score, y=r2_dir, fill = Type, group=P_threshold))
  p <- p +
    geom_bar(stat = "identity", aes(colour = Type), position = "dodge") +
    geom_text(data=subset(top_five_pathways_all, p < 0.05),
              aes(x=score,y=r2_dir,label=p_value_text, hjust=ifelse(sign(r2_dir)>0, 0, 0)), angle = 90, position = position_dodge(width = 1), size = 2.9)
  
  p <- p + scale_x_discrete(labels=c("SCORE_0.05" = "Whole Genome PRS", "SCORE_0.5" = "Whole Genome PRS"))
  p <- p + scale_y_continuous(expand = expand_scale(mult = c(0.2,.6)))
  p <- p + facet_grid(. ~ P_threshold, scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  p <- p + theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1,vjust = 0.5))
  p <- p + ggtitle(top_five_pathways_all$test[1])
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  p <- p + ylab(label = "R2_dir (%)")
  p <- p + xlab(label = "Polygenic risk score")
  
  gglist_r2_dir_pathways[[i]] <- p
  
}




pdf("~/Desktop/ALSPAC_CLOZUK_pathway_whole_genome_comparison.pdf", onefile = TRUE)
#pdf("~/Desktop/ALSPAC_CLOZUK_pathway_whole_genome_comparison_rotated.pdf", onefile = TRUE)
pdf("~/Desktop/ALSPAC_CLOZUK_pathway_whole_genome_comparison_BETAs.pdf", onefile = TRUE)
pdf("~/Desktop/ALSPAC_CLOZUK_pathway_whole_genome_comparison_R2_dir.pdf", onefile = TRUE)

pdf("~/Dropbox/Untitled Folder/ALSPAC_CLOZUK_pathway_whole_genome_genic_comparison_R2_dir_9_pathways.pdf", onefile = TRUE)

#invisible(lapply(gglist_pathways,print))
invisible(lapply(gglist_beta_pathways,print))
invisible(lapply(gglist_r2_dir_pathways,print))
dev.off()

data <- read.csv("~/Downloads/DataJul.csv")

data$rotation[data$Rot.trt %in% c("C2", "S2")]<-"TwoYear"
data$rotation[data$Rot.trt %in% c("C3", "S3", "O3")]<-"ThreeYear"
data$rotation[data$Rot.trt %in% c("C4", "S4", "O4", "A4")]<-"FourYear"

##plot, by rotation #scales = free_x X axis depends on facet
data$rotation <- factor(data$rotation, levels = c("TwoYear", "ThreeYear", "FourYear"))


library(plyr)
data$Rot.Herb.label <- mapvalues(data$Rot.Herb, 
                                 c('C2conv', 'C2low', 'S2conv', 'S2low', 'C3conv', 'C3low',
                                   'O3conv', 'O3low', 'S3conv', 'S3low', 'A4conv', 'A4low',
                                   'C4conv', 'C4low', 'O4conv', 'O4low', 'S4conv', 'S4low'),
                                 c("Corn C", "Corn L", "Soybean C", "Soybean L", 
                                   "Corn C", "Corn L", "Oat C", "Oat L", "Soybean C", 
                                   "Soybean L", "Alfalfa C", "Alfalfa L", "Corn C", "Corn L",
                                   "Oat C", "Oat L", "Soybean C", "Soybean L"))


ggplot(data, aes(Rot.Herb.label, kg.ha, fill=Crop))+
  geom_boxplot()+
  facet_grid(~rotation, scales = "free_x", space="free_x")+
  +
  ggtitle("Weed biomass by plot")+
  theme(plot.title = element_text(size=30, face="bold", vjust=2))+
  xlab("Rotation systems and Herbicide regimes (L = Low herbicide regime, C = Conventional herbicide regime)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Weed dry weight")


### R squared plots ###
# Don't change the code below -- this reformats the data so that it can be inputted to the plot properly

# ensure that the pvalue and r2 columns are classified as numeric by R

results <- ALSPAC_all_pathways
results$p <- as.numeric(results$p)
results$r2_dir <- 100 * (as.numeric(results$r.squared) *
                           (sign(as.numeric(results$estimate))))

### reformat data for plots 
results_p <- results[,c(1,2,6)]
results_r2 <- results[,c(1,2,10)]
results_p_wide <- dcast(results_p, score ~ test , value.var="p")
results_r2_wide <- dcast(results_r2, score ~ test , value.var="r2_dir")

### make plot
### CHANGE PLOTNAME to the name you want for the plot
### the width and height can be adjusted to optimise the aesthetic 
tiff("test.tiff", 
     res = 300, compression = c("lzw"), 
     width = 3000, height = 1500)

# do not change the code below
par(mfrow = c(1,1), "mar" = c(1.5, 3.25, 2.5, 2.5))
mycol = heat.colors(dim(results_r2_wide)[1]);

## npth = number of  p-thresholds
# CHANGE this to the number of p-threshold you used
npth = 2 

### results_r2_wide[,c(2:9)] takes the number of columns with data
### CHANGE : this will vary depending on how many phenotypes of interested you used and would need to be alter 
### should always started at column 2, since the first column are labels for the threshold
### the last column will be 1 plus the number of phenotypes used (here it was 8)
plotmatrix <- as.matrix(results_r2_wide[,c(2:4)])

# do not change the code below
maxp= max(plotmatrix, na.rm = T)
minp= min(plotmatrix, na.rm = T)
minp = min(0, minp)

# CHANGE :  the title of the plot is set by "main = " - change PLOT TITLE to the title that you want 
# CHANGE : ylim sets the limits for the y-axis - depending on the results of your data, thic may need to be altered to make sure everything fits -- currently it is looking at very low R2 values (0.4%) - and allows for bidirectional effects, which might not always be needed if you don't have positive and negative effects (things that increase and decrease risk depending on phenotype).
rip <- barplot(plotmatrix, beside = T, legend = F, col = mycol, 
               ylim = c(-0.4,0.4), 
               main = "TEST_TITLE", 
               cex.main = 1, cex.names = 0.7, 
               axisnames = F, ylab = "",las = 1)

# this may need to be altered to ensure the location of "R-squared (%)" is correct on the plot - this will only be the case if the dimensions of the plot have been alter
# line and padj alter the located of the R-squared
mtext("R-squared (%)", side = 2, cex = 1, las = 1, 
      line = -3, padj = -18.5)

label_names <- colnames (results_r2_wide[,-1])
# CHANGE : the labels are the names of the phenotypes of interest. If multiple phenotypes have been used, make sure you are putting them in the SAME order as they are in the data
# order is determined by the columns of the data frame "results_r2_wide" with the first column being the first print and so on
# CHANGE : the location of the labels on the y-axis is set with "y = " - and this can be changed to have them closer or futher away from the x-axis.   
text(x = colMeans(rip), offset = 1.5, pos = 1, xpd = T , cex = 0.8, srt = 0, 
     y = -0.3,
     labels = label_names)

# CHANGE : The labels below allow for bidirectional effects -- if you don't have bidirection effect just REMOVE these lines.  
# The text can be change to state with it "Increases" or "Decreases" risk, but careful thought must be taken when deciding this effect, keeping in the mind the information from original discovery dataset and the direction of effect within that dataset.  
# location of the text can be altered using padj
mtext("(Positive\nBeta)", 
      padj = -1.5,
      side = 4, line = -1, cex = 0.7, las = 1)
mtext("(Negative\nBeta)", 
      padj = 1.5,
      side = 4, line = -1, cex = 0.7, las = 1)

# do not change this
my.cex = .5
for (cc in 2:9) {
  sig = results_p_wide[,cc];
  sigo = results_p_wide[,cc];
  sig = ifelse(sigo < 0.05, 
               paste("p =", scientific(sigo, digits = 2), sep = " "), "")
  ccl = (cc - 2) * npth
  ripx = rip [(ccl+1):(ccl+npth)]
  sign = ifelse (results_r2_wide[,cc] < 0, 0, results_r2_wide[,cc])
  text(x = ripx, y = sign, sig, srt = 80, cex = 0.7, adj = c(-0.1,0))
}

# CHANGE : the legend states which pvalue thresholds were used - alter this to correctly reflect the thresholds that were used.   
legend("top", title = "P-value threshold in training data", 
       legend = c("p<0.05",
                  "p<0.5"),  
       fill = mycol, ncol = length(results_r2_wide[,1]), cex =.7)

dev.off()


library(plotly)
library(dplyr)

p <- ggplot2::diamonds %>% count(cut, clarity) %>%
  plot_ly(x = ~cut, y = ~n, color = ~clarity)


# Plotly plots 
# run with all brain regions, assigning each to a huge table eg: with i = 1 and 2
testing_pathways <- rbind(testing_1_pathways, top_five_pathways_all)
testing_pathways$Typetwo <- c("Pathway1", "Pathway2", "Pathway3", "Pathway4", "Pathway5", "Whole_genome","Pathway1", "Pathway2", "Pathway3", "Pathway4", "Pathway5", "Whole_genome", "Pathway1", "Pathway2", "Pathway3", "Pathway4", "Pathway5", "Whole_genome", "Pathway1", "Pathway2", "Pathway3", "Pathway4", "Pathway5", "Whole_genome")
testing_pathway_1 <- testing_pathways[grep("Pathway1", testing_pathways$Typetwo),]
Pathway_1_rsquared <- as.numeric(testing_pathway_1$r.squared)
test <- testing_pathway_1$test
testing_pathway_2 <- testing_pathways[grep("Pathway2", testing_pathways$Typetwo),]
Pathway_2_rsquared <- as.numeric(testing_pathway_2$r.squared)
Pathway_1_pathways <- as.character(testing_pathway_1$score)
Pathway_2_pathways <- as.character(testing_pathway_2$score)
# Need the p_value as well, can add in later
# Make two separate graphs, one for each p_value threshold
final_table <- data.frame(test,Pathway_1_rsquared,Pathway_1_pathways, Pathway_2_rsquared, Pathway_2_pathways)
final_table_test <- final_table[c(1,3),] # seperate into separate p value thresholds

p <- plot_ly(final_table_test, x= ~test, y= ~Pathway_1_rsquared, type = 'bar', name = "Pathway one", 
             text = ~Pathway_1_pathways, textposition = 'auto', marker = list(color = 'rgb(49,130,189)')) %>% 
  add_trace(y= ~Pathway_2_rsquared, name = 'Pathway two', text = ~Pathway_2_pathways, textposition = 'auto', marker = list(color = 'rgb(49,125,189)')) %>% 
  layout(xaxis = list(title = "", tickangle = -45), yaxis = list(title = ""), margin = list(b=200), barmode = 'group')
