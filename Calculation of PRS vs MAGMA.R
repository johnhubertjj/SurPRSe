library(zoom)
library(ggplot2)
library(data.table)

setwd("/Volumes/NO NAME/LOG10 graphs/")
wd <- getwd()
sig <- c("0.0001","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")
sig2 <- c(1e-4,1e-3,1e-2,0.05,0.1,0.2,0.3,0.4,0.5)
correlations <- rep(0,length(sig2))

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

gglist <- list()

for (i in 1:length(sig)) {
MAGMA <- fread(paste0(sig[i], "gene_annotation_for_CLOZUK_test.genes.out"))
PRS <- fread(paste0(sig2[i], "CLOZUK_PGC_PRS_residuals_with_genes_using_fam_file.txt"))
names(PRS) <- c("P","GENE")

Merged_table <- merge(PRS,MAGMA,by = 'GENE')

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

plot.subtitle <- paste0("correlation = ", format(round(correlations[i], 2), nsmall = 2))

gglist[[i]] <- ggplot(Merged_table, aes(graph_PRS,graph_magma)) + geom_point() + geom_abline() + theme(axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
gglist[[i]] <- gglist[[i]] + xlim(0,8) + ylim(0,8) + labs (x = expression(PRS~~-log[10](italic(p))), y = expression(MAGMA~~-log[10](italic(p)))) + ggtitle(sig[i], subtitle = plot.subtitle)
gglist[[i]] <- gglist[[i]] + theme(plot.title = element_text(family = "", color="#666666", face="bold", size=32, hjust=0))

finding_out_magma_SNPS <- subset(Merged_table, (logp.y > 2))
finding_out_the_rest <- subset(Merged_table, (logp.y <= 2))

sd_interesting <- sd(finding_out_magma_SNPS$NPARAM)
sd_compare <- sd(finding_out_the_rest$NPARAM)
mean_interesting <- mean(finding_out_magma_SNPS$NPARAM)
mean_compare <- mean(finding_out_the_rest$NPARAM)

SNPs_per_genes_MAGMA_sd_interesting[i] <- sd_interesting
SNPs_per_genes_MAGMA_mean_interesting[i] <- mean_interesting
SNPs_per_genes_MAGMA_mean_compare[i] <- mean_compare
SNPs_per_genes_MAGMA_sd_compare[i] <- sd_compare
number_of_SNPs_MAGMA[i] <- sum(MAGMA$NPARAM)
}

gglist[[i]]
mypath <- paste0(wd,sig[i],"PRSvsMAGMA_CLOZUK.png")
ggsave(filename = mypath,plot = plots)
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

