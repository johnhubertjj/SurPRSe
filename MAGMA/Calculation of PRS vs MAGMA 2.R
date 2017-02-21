library(zoom)
library(ggplot2)
setwd("/Volumes/NO NAME/LOG10 graphs/")
wd <- getwd()
sig <- c("0.0001","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")
sig2 <- c(1e-4,1e-3,1e-2,0.05,0.1,0.2,0.3,0.4,0.5)
correlations <- rep(0,length(sig2))
par(mfrow=c(3,3))

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

par(mfrow=c(3,3))

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

plots <- ggplot(Merged_table, aes(graph_PRS,graph_magma)) + geom_point() + geom_abline() + theme(axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
plots <- plots + xlim(0,8) + ylim(0,8) + labs (x = expression(PRS~~-log[10](italic(p))), y = expression(MAGMA~~-log[10](italic(p)))) + ggtitle(paste0 (sig[i], "correlations = ", correlations[i]))
plots                                                                                                                                
}

ggd.qqplot(Merged_table$P.y,Merged_table$P.x,title)
ggd.qqplot.limit(Merged_table$P.y,Merged_table$P.x,title)
ggd.qqplot.in.R(Merged_table$P.y,Merged_table$P.x,title)
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


