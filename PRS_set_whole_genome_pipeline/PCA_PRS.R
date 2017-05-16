
data(iris)
library(data.table)

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


Schizophrenia <- fread("E:/PGC_CLOZUK_output/PRS_scoring/PGC_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")
Bipolar <- fread("E:/BIP_CLOZUK_output/PRS_scoring/BIP_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")
Educational_attainment <- fread("E:/EDU_main_CLOZUK_output/PRS_scoring/EDU_main_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")
PGC_MDD <- fread("E:/PGC_MDD_CLOZUK_output/PRS_scoring/PGC_MDD_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")
BIPvsSCZ <- fread("E:/BIPvsSCZ_CLOZUK_output/PRS_scoring/BIPvsSCZ_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")
Neuroticism <- fread("E:/Neurot_Assoc_Biobank_CLOZUK_output/PRS_scoring/Neurot_Assoc_Biobank_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")

covariates <- fread("F:/Stationary_data/CLOZUK2.r7.select2PC.eigenvec.txt")
fam2 <- fread("F:/Stationary_data/CLOZUK.r7.GWAS_IDs.fam")
colnames(fam2) <- c("FID","IID","PID","MID","Sex","PHENO")

PRS.profiles.1 <- Neuroticism

PRS.Profiles.with.covariates <- merge(covariates,PRS.profiles.1,by.x="FID", by.y="FID", all = F)
PRS.Profiles.with.covariates <- merge(PRS.Profiles.with.covariates, fam2, by.x = "FID", by.y = "FID", all = F)

#  res$model[i]<-sig[i]

# Calculate model including covariates against the polygenic risk score
model0<-glm(SCORE~PC1+PC2+PC3+PC4+PC5+PC6+PC9+PC11+PC12+PC13+PC19, family = gaussian, data = PRS.Profiles.with.covariates)
m1<-mean(residuals(model0))
sd1<-sd(residuals(model0))

# Calculate the Normalised score
PRS.Profiles.with.covariates$NORMSCORE<-(residuals(model0)-m1)/sd1

#             hist(PRS.profiles$NORMSCORE)

# change the Phenotypes so that they will work in a binary model
PRS.Profiles.with.covariates$PHENO.y <- PRS.Profiles.with.covariates$PHENO.y - 1

merge1 <- PRS.Profiles.with.covariates
merge2 <- PRS.Profiles.with.covariates
merge3 <- PRS.Profiles.with.covariates
merge4 <- PRS.Profiles.with.covariates
merge5 <- PRS.Profiles.with.covariates
merge6 <- PRS.Profiles.with.covariates

PCA_matrix <- fread("/Volumes/PhD_storage/PRS_cross_disorder_table.csv")

Groups_of_individuals <- c("CLOZUK","COGS","CRESTAR1", "CRESTAR2", "CRESTAR3", "1958BC", "BLOOD", "GERAD", "CON_GS", "HYWEL","POBI","QIMR","T1DGC","TEDS", "TWINSUK","WTCCC")
PCA_matrix$Colours <- "NA"

for (i in 1:length(Groups_of_individuals)){
  Current_integers <-  PCA_matrix[,.I[grep(Groups_of_individuals[i], Individuals) ]]
  PCA_matrix$Colours[Current_integers] <- Groups_of_individuals[i]
} 


#PCA_matrix <- data.frame(data = merge1$FID,merge1$PHENO.y, merge1$NORMSCORE, merge2$NORMSCORE, merge3$NORMSCORE, merge4$NORMSCORE, merge5$NORMSCORE, merge6$NORMSCORE)
library(FactoMineR)

#names(PCA_matrix) <- c("Individuals","PHENOTYPE", "Schizophrenia", "Bipolar", "Educational_attainment", "PGC_MDD", "BIPvsSCZ", "Neuroticism")

#write.csv(PCA_matrix, file = "E:/PRS_cross_disorder_table.csv", col.names = T, row.names = F)

PCA_matrix_df <- as.data.frame(PCA_matrix)
PCA_matrix_df_cases <- as.data.frame(PCA_matrix[PHENOTYPE==1])

testing <- prcomp(PCA_matrix_df[3:8],center = T, scale. = T)
testing2 <- prcomp(PCA_matrix_df_cases[3:8],center = T, scale. = T)

rawLoadings_cases <- testing2$rotation[,1:6] %*% diag(testing2$sdev, 6, 6)
rotatedLoadings <- varimax(rawLoadings_cases)$loadings
invloadings <- t(pracma::pinv(rotatedLoadings))
scores <- scale(PCA_matrix_df_cases[,3:8]) %*% invloadings

colnames(scores) <-c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6")
column_names <- colnames(scores)

varimax_list <- list()
for ( i in 1:15){
  a <- scores[,plots_index[1,i]]
  b <- scores[,plots_index[2,i]]
  
  p <- qplot(a, b, main = "PCA on PRS Cases with Varimax", xlab = column_names[plots_index[1,i]], ylab = column_names[plots_index[2,i]])
  
  varimax_list[[i]] <- p
}


multiplot(plotlist = varimax_list, cols = 3)

biplot(testing)

a <- cov(PCA_matrix_df[3:8])
eigenPRS <- eigen(a)


library(devtools)
install_github("ggbiplot", "vqv")

plots_index <- combn(1:6,2)
e <- new.env()
library(ggbiplot)

gglist <- list()
for ( i in 1:15){

g <- ggbiplot(testing2, obs.scale = 1, var.scale = 1, ellipse = F, choices = c(plots_index[1,i],plots_index[2,i]),
              circle = F ,alpha = 1)
g <- g + scale_color_discrete(name="Phenotype_CLOZUK.BGE")
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

gglist[[i]] <- g
}

multiplot(plotlist = gglist, cols = 3)
screeplot(testing2, main="Scree Plot", xlab="Components",ylim = c(0,1.4))


autoplot(prcomp(new_table2), data = new_table2)
d.factanal <- factanal(testing, factors = 6,rotation = "varimax")

plot(d.factanal$loadings[,1:2], typ="n")
text(d.factanal$loadings[,1:2],labels=names(new.data.frame[,2:6]),cex=.7) # add variable names

library(ggfortify)
library(cluster)
autoplot(prcomp(new.data.frame))

library(corrplot)
corrplot(cor(testing[,1:100]), order = "hclust", tl.col='black', tl.cex=.75) 

plot(PCA_individ$rotation[,1], PCA_individ$rotation[,2])

new.data.frame.cases <- new.data.frame[new.data.frame$pheno == 1,]

ncomp <- 5
pca_iris_rotated <- psych::principal(new.data.frame.cases[,2:6], rotate="varimax", nfactors=ncomp, scores=TRUE)
print(pca_iris_rotated$scores)

new_table$colour <- "black"
new_table$colour[new_table$PHENO.y==1]="red"
testing_lol <-pca_iris_rotated$scores
plot(pca_iris_rotated$scores[,4], pca_iris_rotated$scores[,1]),col = new_table$colour)
plot(pca_iris_rotated$scores[,4], pca_iris_rotated$scores[,2])

testingy <- princomp(new.data.frame.cases[2:7])


library(FactoMineR)
data(decathlon)
res.pca = PCA(new.data.frame[,2:6], scale.unit=TRUE, ncp=5, graph=T)

d_stan = as.data.frame(scale(new.data.frame[,2:6]))
res1b = factanal(d_stan, factors = 2, rotation = "none", na.action = na.omit)
plot(res1b$loadings)

### Plot loadings against one another
load = res1b$loadings[,1:2]
plot(load, type="n") # set up plot 
text(load,labels=names(d_stan),cex=.7) # add variable names

irisX <- iris[1:nrow(iris), 1:4]      # Iris data
ncomp <- 2

pca_iris_rotated <- psych::principal(irisX, rotate="none", nfactors=ncomp, scores=TRUE)
print(pca_iris_rotated$scores[1:5,])  # Scores returned by principal()
plot(pca_iris_rotated$scores)

PCA_matrix_cases <- PCA_matrix[PCA_matrix$PHENOTYPE == 1,]
PCA_matrix_controls <- PCA_matrix[PCA_matrix$PHENOTYPE == 0,]


my.prc <- prcomp(PCA_matrix[,-1:-2], scale. = T, center = T)
my.prc_cases <- prcomp(PCA_matrix_cases[,-1:-2], scale. = F, center = T) 
my.prc_controls<- prcomp(PCA_matrix_controls[,-1:-2], scale. = T, center = T) 

rawLoadings_cases <- my.prc_cases$rotation[,3:8] %*% diag(my.prc_cases$sdev, 6, 6)
rotatedLoadings <- varimax(rawLoadings_cases)$loadings
invloadings <- t(pracma::pinv(rotatedLoadings))
scores <- scale(PCA_matrix_cases[,3:8]) %*% invloadings


iris

my.abs     <- abs(cor(PCA_matrix_cases[,-1:-2]))
my.colors  <- dmat.color(my.abs)
my.ordered <- order.single(cor(PCA_matrix_cases[,-1:-2]))
cpairs(PCA_matrix_cases[,-1:-2], my.ordered, panel.colors=my.colors, gap=0.5)

my.abs     <- abs(cor(PCA_matrix_controls[,-1:-2]))
my.colors  <- dmat.color(my.abs)
my.ordered <- order.single(cor(PCA_matrix_controls[,-1:-2]))
cpairs(PCA_matrix_controls[,-1:-2], my.ordered, panel.colors=my.colors, gap=0.5)

library(lattice)
library(mclust)

dat.em <- mclustBIC(PCA_matrix_df[,c(3:8)]) 
splom(as.data.frame(dat.pca$x), 
      col=summary(dat.em,data=dat)$classification, cex=2,pch='*')


