ibrary(zoom)
library(ggplot2)
library(data.table)
library(grid)
library(gridExtra)

setwd("/Volumes/NO NAME/LOG10 graphs/")
wd <- getwd()
sig <- c("0.0001","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")
sig2 <- c(1e-4,1e-3,1e-2,0.05,0.1,0.2,0.3,0.4,0.5)
correlations <- rep(0,length(sig2))

# environment
e <- new.env()
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## qqplot function with limits
ggd.qqplot.limit <- function(pvector1, pvector2, main=NULL, ...) {
  o = -log10(pvector1)
  e = -log10(pvector2)
  mypath <- paste0(wd,sig[i],"PRSvsMAGMA_CLOZUK.png")
  png(file = mypath)
  plot(e,o,pch=1,cex=1, main=main, ...,
       xlab=expression(PRS~~-log[10](italic(p))),
       ylab=expression(MAGMA~~-log[10](italic(p))),
       xlim=c(0,20), ylim=c(0,20))
  lines(e,e,col="red")
  dev.off()
}
## qqfunction without limits
ggd.qqplot <- function(pvector1, pvector2, main=NULL, ...) {
  o = -log10(pvector1)
  e = -log10(pvector2)
  mypath <- paste0(wd,"/",sig[i],"PRSvsMAGMA_CLOZUK.png")
  png(file = mypath)
  plot(e,o,pch=1,cex=1, main=main, ...,
       xlab=expression(PRS~~-log[10](italic(p))),
       ylab=expression(MAGMA~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
  dev.off()
}

## qqfuncton in R environ.
ggd.qqplot.in.R <- function(pvector1, pvector2, main=NULL, ...) {
  o = -log10(pvector1)
  e = -log10(pvector2)
  plot(e,o,pch=1,cex=1, main=main, ...,
       xlab=expression(PRS~~-log[10](italic(p))),
       ylab=expression(MAGMA~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  ?lines(e,e,col="red")
}

number_of_combined_genes <- rep(0,length(sig))
number_genes_PRS <-rep(0,length(sig))
number_genes_MAGMA <- rep(0,length(sig))
SNPs_per_genes_MAGMA_sd_interesting <- rep(0, length(sig))
SNPs_per_genes_MAGMA_mean_interesting <- rep(0, length(sig))
SNPs_per_genes_MAGMA_sd_compare <- rep(0, length(sig))
SNPs_per_genes_MAGMA_mean_compare <- rep(0, length(sig))
number_of_SNPs_MAGMA <- rep(0, length(sig))

SNPs_per_genes_MAGMA_sd_interesting2 <- rep(0, length(sig))
SNPs_per_genes_MAGMA_mean_interesting2 <- rep(0, length(sig))
SNPs_per_genes_MAGMA_sd_compare2 <- rep(0, length(sig))
SNPs_per_genes_MAGMA_mean_compare2 <- rep(0, length(sig))

gglist <- list()

testing_PRS_VS_MAGMA <- function() {
for (i in 1:length(sig)) {
  MAGMA <- fread(paste0("./output/MAGMA_set_analysis/",sig[i], "gene_annotation_for_CLOZUK_without_selecting_genes_magma.genes.out"))
   #PRS <- fread(paste0("~/Dropbox/PhDwork-PRS/",sig2[i], "CLOZUK_PGC_PRS_residuals_with_genes_using_fam_file_whole_genome_more_than_one_SNP.txt"))
   PRS <- fread(paste0("~/Dropbox/PhDwork-PRS/",sig2[i], "CLOZUK_PGC_PRS_residuals_with_genes_using_fam_file_whole_genome.txt"))
   names(PRS) <- c("P","GENE")
   #names(PRS) <- c("GENE", "P")
  
  Merged_table <- merge(PRS,MAGMA,by = 'GENE')
  Merged_table <- Merged_table[NSNPS > 1]
  
  Merged_table[,logp.x := -log10(P.x)]
  Merged_table[,logp.y := -log10(P.y)]
  
  title <- sig[i]
  PRS.pvalue <- Merged_table$P.x
  MAGMA.P <- Merged_table$P.y
  
  correlations[i] <- cor(PRS.pvalue, MAGMA.P)
  number_of_combined_genes[i] <- nrow(Merged_table)
  number_genes_PRS[i] <- length(which(!is.na(PRS$GENE)))
  number_genes_MAGMA[i] <- nrow(MAGMA)
  
  graph_PRS <- -log10(PRS.pvalue)
  graph_magma <- -log10(MAGMA.P)
  
  correlations[i] <- signif(correlations[i], digits = 2)
  plot.subtitle <- paste0("correlation = ", as.character(correlations[i]))
  
  p1 <- ggplot(Merged_table, aes_string(x = "logp.x", y = "logp.y"))+ geom_point() + geom_abline() + theme(axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  p1 <- p1 + xlim(0,8) + ylim(0,8) + labs (x = expression(PRS~~-log[10](italic(p))), y = expression(MAGMA~~-log[10](italic(p)))) + ggtitle(sig[i], subtitle = plot.subtitle)
  p1 <- p1 + theme(plot.title = element_text(family = "", color="#666666", face="bold", size=32, hjust=0))
  
  gglist[[i]] <- p1
  assign(paste0("plot",sig[i]), p1, envir = e)
  
  
  finding_out_magma_SNPS <- subset(Merged_table, (logp.y > 10 & logp.x < 10))
  finding_out_the_rest <- subset(Merged_table, (!(logp.y > 10 & logp.x < 10)))
  #browser()
  sd_interesting <- sd(finding_out_magma_SNPS$NPARAM)
  sd_compare <- sd(finding_out_the_rest$NPARAM)
  mean_interesting <- mean(finding_out_magma_SNPS$NPARAM)
  mean_compare <- mean(finding_out_the_rest$NPARAM)
  
  sd_interesting2 <- sd(finding_out_magma_SNPS$NSNPS)
  sd_compare2 <- sd(finding_out_the_rest$NSNPS)
  mean_interesting2 <- mean(finding_out_magma_SNPS$NSNPS)
  mean_compare2 <- mean(finding_out_the_rest$NSNPS)
  
  SNPs_per_genes_MAGMA_sd_interesting2[i] <- sd_interesting2
  SNPs_per_genes_MAGMA_mean_interesting2[i] <- mean_interesting2
  SNPs_per_genes_MAGMA_mean_compare2[i] <- mean_compare2
  SNPs_per_genes_MAGMA_sd_compare2[i] <- sd_compare2
  
  SNPs_per_genes_MAGMA_sd_interesting[i] <- sd_interesting
  SNPs_per_genes_MAGMA_mean_interesting[i] <- mean_interesting
  SNPs_per_genes_MAGMA_mean_compare[i] <- mean_compare
  SNPs_per_genes_MAGMA_sd_compare[i] <- sd_compare
  number_of_SNPs_MAGMA[i] <- sum(MAGMA$NSNPS)
  #browser()
  #print(gglist[[i]])
  
  #  mypath <- paste0(wd,sig[i],"PRSvsMAGMA_CLOZUK_without_one_SNP.png")
  #  mypath <- paste0(wd,sig[i],"PRSvsMAGMA_CLOZUK_no_boundaries.png")
  #  mypath <- paste0(wd,sig[i],"PRSvsMAGMA_CLOZUK.png")
  # ggsave(filename = mypath,plot = gglist[[i]])
}
assign("SNPs_per_genes_MAGMA_sd_interesting2",SNPs_per_genes_MAGMA_sd_interesting2, envir = e)  
assign("SNPs_per_genes_MAGMA_sd_compare2",SNPs_per_genes_MAGMA_sd_compare2, envir = e)
assign("SNPs_per_genes_MAGMA_sd_compare",SNPs_per_genes_MAGMA_sd_compare, envir = e)
assign("SNPs_per_genes_MAGMA_sd_interesting",SNPs_per_genes_MAGMA_sd_interesting, envir = e)
assign("SNPs_per_genes_MAGMA_mean_compare2",SNPs_per_genes_MAGMA_mean_compare2, envir = e)
assign("SNPs_per_genes_MAGMA_mean_compare",SNPs_per_genes_MAGMA_mean_compare, envir = e)
assign("SNPs_per_genes_MAGMA_mean_interesting",SNPs_per_genes_MAGMA_mean_interesting, envir = e)
assign("SNPs_per_genes_MAGMA_mean_interesting2",SNPs_per_genes_MAGMA_mean_interesting2, envir = e)

names_for_matrix <- list()
names_for_matrix[[2]] <- sig
log_names <- c("Outlier_data(logp more than 10)","Correlated data(logp less than 10)")
names_for_matrix[[1]] <- log_names

Mean_for_NSNPS <- matrix(c(e$SNPs_per_genes_MAGMA_mean_interesting2,SNPs_per_genes_MAGMA_mean_compare2), nrow = 2, byrow = T,dimnames = names_for_matrix)
SD_for_NSNPS <- matrix(c(e$SNPs_per_genes_MAGMA_sd_interesting2,SNPs_per_genes_MAGMA_sd_compare2), nrow = 2, byrow = T,dimnames = names_for_matrix)
Mean_for_NPARAMS <- matrix(c(e$SNPs_per_genes_MAGMA_mean_interesting,SNPs_per_genes_MAGMA_mean_compare), nrow = 2, byrow = T, dimnames = names_for_matrix)
SD_for_NPARAMS <- matrix(c(e$SNPs_per_genes_MAGMA_sd_interesting,SNPs_per_genes_MAGMA_sd_compare), nrow = 2, byrow = T, dimnames = names_for_matrix)
test_list <- list(Mean_for_NSNPS,SD_for_NSNPS,Mean_for_NPARAMS,SD_for_NPARAMS)
names_of_tables <- c("Mean_for_NSNPS","SD_for_NSNPS","Mean_for_NPARAMS","SD_for_NPARAMS")
for (i in 1:4){
  write.csv(test_list[[i]], file = paste0(names_of_tables[i],".csv"))
} 

assign("Mean_for_NSNPs", Mean_for_NSNPS, envir = .GlobalEnv)
assign("SD_for_NSNPs", SD_for_NSNPS, envir = .GlobalEnv)
assign("Mean_for_NPARAMS", Mean_for_NPARAMS, envir = .GlobalEnv)
assign("SD_for_NPARAMS", SD_for_NPARAMS, envir = .GlobalEnv)

layout <- matrix(c(1, 2, 3, 4, 5, 6, 7,8,9), nrow = 3, byrow = TRUE)
multiplot(plotlist = gglist, cols = 3, layout = layout)
}

layout <- matrix(c(1, 2, 3, 4, 5, 6, 7,8,9), nrow = 3, byrow = TRUE)

for (i in 1:9) {
  mypath <- paste0(wd,sig[i],"PRSvsMAGMA_CLOZUK_without_one_SNP.png")
  ggsave(filename = mypath,plot = gglist[[i]])
}

multiplot(gglist,col =3)
mypath <- paste0(wd,sig[i],"PRSvsMAGMA_CLOZUK.png")
png(file = mypath)
plots
dev.off()
assign(paste0("plot", i), plots)
}

multiplot(gglist, cols = 3)

ggd.qqplot(Merged_table$P.y,Merged_table$P.x,title)
ggd.qqplot.limit(Merged_table$P.y,Merged_table$P.x,title)
ggd.qqplot.in.R(Merged_table$P.y,Merged_table$P.x,title)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# http://GettingGeneticsDone.blogspot.com/
# See http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Define the function


## testing with ggplot

ggplot(Merged_table, aes(graph_PRS,graph_magma)) + geom_point() + geom_abline()

par(mfrow=c(3,3))

for (i in 1:9){
  plots <- ggplot(Merged_table, aes(graph_PRS,graph_magma)) + geom_point() + geom_abline() + theme(axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  plots <- plots + xlim(0,8) + ylim(0,8) + labs (x = expression(PRS~~-log[10](italic(p))), y = expression(MAGMA~~-log[10](italic(p)))) + ggtitle(paste0 (sig[i], "correlations = ", correlations[i])
                                                                                                                                                 plots
}

plots <- plots + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))


# plots <- plots + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
# ggplot(Merged_table, aes(graph_PRS,graph_magma)) + geom_point() + geom_abline() + theme(axis.line.x = element_line(colour = "black", size = 3))
ggtitle(bquote(atop(.(sig[i]), atop(italic(.(plot.subtitle)), ""))))


hist(finding_out_the_rest$NSNPS, col=rgb(0,0,1,0.5),breaks = 100)
hist(finding_out_magma_SNPS$NSNPS, col=rgb(1,0,0,0.5), breaks = 100)