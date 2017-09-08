
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


#Schizophrenia <- fread("/Volumes/PhD_storage/PGC_CLOZUK_output/PRS_scoring/PGC_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")
Schizophrenia <- fread("~/Dropbox/PRS_scores_different_summary_stats_datasets/PGC_CLOZUK_whole_genome_significance_threshold_at_0.1.profile")
#Bipolar <- fread("/Volumes/PhD_storage/BIP_CLOZUK_output/PRS_scoring/BIP_CLOZUK_whole_genome_significance_threshold_at_0.4.profile")
Bipolar <- fread("~/Dropbox/PRS_scores_different_summary_stats_datasets/BIP_CLOZUK_whole_genome_significance_threshold_at_0.4.profile")
#Educational_attainment <- fread("/Volumes/PhD_storage/EDU_main_CLOZUK_output/PRS_scoring/EDU_main_CLOZUK_whole_genome_significance_threshold_at_0.15.profile")
Educational_attainment <- fread("~/Dropbox/PRS_scores_different_summary_stats_datasets/EDU_main_CLOZUK_whole_genome_significance_threshold_at_0.15.profile")
#PGC_MDD <- fread("/Volumes/PhD_storage/PGC_MDD_CLOZUK_output/PRS_scoring/PGC_MDD_CLOZUK_whole_genome_significance_threshold_at_0.45.profile")
PGC_MDD <- fread("~/Dropbox/PRS_scores_different_summary_stats_datasets/PGC_MDD_CLOZUK_whole_genome_significance_threshold_at_0.45.profile")
#BIPvsSCZ <- fread("/Volumes/PhD_storage/BIPvsSCZ_CLOZUK_output/PRS_scoring/BIPvsSCZ_CLOZUK_whole_genome_significance_threshold_at_0.35.profile")
BIPvsSCZ <- fread("~/Dropbox/PRS_scores_different_summary_stats_datasets/BIPvsSCZ_CLOZUK_whole_genome_significance_threshold_at_0.35.profile")
#Neuroticism <- fread("/Volumes/PhD_storage/Neurot_Assoc_Biobank_CLOZUK_output/PRS_scoring/Neurot_Assoc_Biobank_CLOZUK_whole_genome_significance_threshold_at_0.15.profile")
Neuroticism <- fread("~/Dropbox/PRS_scores_different_summary_stats_datasets/Neurot_Assoc_Biobank_CLOZUK_whole_genome_significance_threshold_at_0.15.profile")
IQ <- fread("~/Dropbox/IQ_2017_PROFILES/PRS_scoring/IQ_GWAS_2017_CLOZUK_whole_genome_significance_threshold_at_0.2.profile")

Neuropsychiatric_datasets <- list(Schizophrenia,Bipolar,Educational_attainment, PGC_MDD, BIPvsSCZ, Neuroticism, IQ)

#Neuropsychiatric_datasets <- list(Schizophrenia,Bipolar,Educational_attainment, PGC_MDD, BIPvsSCZ, Neuroticism)

covariates <- fread("~/Dropbox/Stationary_data/CLOZUK2.r7.select2PC.eigenvec.txt")
fam2 <- fread("~/Dropbox/Stationary_data/CLOZUK.r7.GWAS_IDs.fam")
colnames(fam2) <- c("FID","IID","PID","MID","Sex","PHENO")

Groups_to_keep <- c("CLOZUK","COGS","CRESTAR1", "CRESTAR2", "CRESTAR3","POBI","T1DGC")
#Groups_to_keep <- c("CLOZUK","COGS","CRESTAR1", "CRESTAR2", "CRESTAR3","TWINSUK","1958BC")

for (i in 1:7){
setkey(Neuropsychiatric_datasets[[i]],FID)
 #Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
                      #Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
 


PRS.profiles.1 <- Neuropsychiatric_datasets[[i]]

PRS.Profiles.with.covariates <- merge(covariates,PRS.profiles.1, by.x="FID", by.y="FID", all = F)
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

Neuropsychiatric_datasets[[i]] <- PRS.Profiles.with.covariates

#Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
#Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
}

for (i in 1:6){
  setkey(Neuropsychiatric_datasets[[i]],FID)
  # Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
  #                      Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
  
  
  
  PRS.profiles.1 <- Neuropsychiatric_datasets[[i]]
  
  PRS.Profiles.with.covariates <- merge(covariates,PRS.profiles.1, by.x="FID", by.y="FID", all = F)
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
  
  Neuropsychiatric_datasets[[i]] <- PRS.Profiles.with.covariates
  
  Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
                                                                        Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]
}

PCA_matrix_2 <- data_frame(Neuropsychiatric_datasets[[1]]$FID,Neuropsychiatric_datasets[[1]]$PHENO.y,Neuropsychiatric_datasets[[1]]$NORMSCORE, Neuropsychiatric_datasets[[2]]$NORMSCORE,Neuropsychiatric_datasets[[3]]$NORMSCORE, Neuropsychiatric_datasets[[4]]$NORMSCORE,Neuropsychiatric_datasets[[5]]$NORMSCORE, Neuropsychiatric_datasets[[6]]$NORMSCORE, Neuropsychiatric_datasets[[7]]$NORMSCORE)
PCA_matrix_2 <- data_frame(Neuropsychiatric_datasets[[1]]$FID,Neuropsychiatric_datasets[[1]]$PHENO.y,Neuropsychiatric_datasets[[1]]$NORMSCORE, Neuropsychiatric_datasets[[2]]$NORMSCORE,Neuropsychiatric_datasets[[3]]$NORMSCORE, Neuropsychiatric_datasets[[4]]$NORMSCORE,Neuropsychiatric_datasets[[5]]$NORMSCORE, Neuropsychiatric_datasets[[6]]$NORMSCORE)

test <- colnames(PCA_matrix)
colnames(PCA_matrix_2) <- test

merge1 <- PRS.Profiles.with.covariates
merge2 <- PRS.Profiles.with.covariates
merge3 <- PRS.Profiles.with.covariates
merge4 <- PRS.Profiles.with.covariates
merge5 <- PRS.Profiles.with.covariates
merge6 <- PRS.Profiles.with.covariates

PCA_matrix <- fread("/Volumes/PhD_storage/PRS_cross_disorder_table_optimised_thresholds.csv")
#PCA_matrix <- fread("E:/PRS_cross_disorder_table_optimised_thresholds.csv")

Groups_of_individuals <- c("CLOZUK","COGS","CRESTAR1", "CRESTAR2", "CRESTAR3", "1958BC", "BLOOD", "GERAD", "CON_GS", "HYWEL","POBI","QIMR","T1DGC","TEDS","TWINSUK","WTCCC")

PCA_matrix$Colours <- "NA"

for (i in 1:length(Groups_of_individuals)){
  Current_integers <-  PCA_matrix[,.I[grep(Groups_of_individuals[i], PCA_matrix$Individuals) ]]
  PCA_matrix$Colours[Current_integers] <- Groups_of_individuals[i]
} 


#PCA_matrix <- data.frame(data = merge1$FID,merge1$PHENO.y, merge1$NORMSCORE, merge2$NORMSCORE, merge3$NORMSCORE, merge4$NORMSCORE, merge5$NORMSCORE, merge6$NORMSCORE)
library(FactoMineR)


names(PCA_matrix_2) <- c("Individuals","PHENOTYPE", "Schizophrenia", "Bipolar", "Educational_attainment", "PGC_MDD", "BIPvsSCZ", "Neuroticism", "IQ")

#write.csv(PCA_matrix, file = "/Volumes/PhD_storage/PRS_cross_disorder_table_optimised_thresholds.csv", col.names = T, row.names = F)

PCA_matrix_df <- as.data.frame(PCA_matrix)
PCA_matrix_df_cases <- as.data.frame(PCA_matrix[PHENOTYPE==1])

testing <- prcomp(PCA_matrix_df[3:8],center = T)
testing_reduce_controls <- prcomp(PCA_matrix_2[3:8],center = T)
testing_reduce_controls <- prcomp(PCA_matrix_2[3:9],center = T)

testing_wo_SCZ <- prcomp(PCA_matrix_df[4:8],center = T)
testing_wo_BIP <- prcomp(PCA_matrix_df[c(3,5:8)],center = T)
testing_wo_EDU <- prcomp(PCA_matrix_df[c(3:4,6:8)],center = T)
testing_wo_MDD <- prcomp(PCA_matrix_df[c(3:5,7:8)],center = T)
testing_wo_BIPvsSCZ <- prcomp(PCA_matrix_df[c(3:6,8)],center = T)
testing_wo_Neurot <- prcomp(PCA_matrix_df[c(3:7)],center = T)

testing_princomp <- princomp(PCA_matrix_df[3:8],center = T)
testing_experiments <- prcomp(PCA_matrix_df_cases[c(3:7)],center = T)
testing2 <- prcomp(PCA_matrix_df_cases[3:8], center = T)

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

for (i in 1:15){

g <- ggbiplot(testing, obs.scale = 1, var.scale = 1, ellipse = F, choices = c(plots_index[1,i],plots_index[2,i]),
              circle = F ,groups = as.factor(PCA_matrix_df$PHENOTYPE), alpha = 1)
g <- g + scale_color_discrete(name="Phenotype_CLOZUK.BGE")
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

gglist[[i]] <- g
}

# Testing_with_different_controls

library(devtools)
install_github("ggbiplot", "vqv")

plots_index <- combn(1:7,2)
e <- new.env()
library(ggbiplot)

gglist <- list()

for (i in 1:20){
  
  g <- ggbiplot(testing_reduce_controls, obs.scale = 1, var.scale = 1, ellipse = F, choices = c(plots_index[1,i],plots_index[2,i]),
                circle = F ,groups = as.factor(PCA_matrix_2$PHENOTYPE), alpha = 1)
  g <- g + scale_color_discrete(name="Phenotype_CLOZUK.BGE")
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
  
  gglist[[i]] <- g
}

multiplot(plotlist = gglist, cols = 4)
screeplot(testing2, main="Scree Plot", xlab="Components",ylim = c(0,1.4))


plots_index <- combn(1:5,2)

gglist <- list()
for ( i in 1:10){
  
  g <- ggbiplot(testing_experiments, obs.scale = 1, var.scale = 1, ellipse = F, choices = c(plots_index[1,i],plots_index[2,i]),
                circle = F ,groups = as.factor(PCA_matrix_df$PHENOTYPE), alpha = 1)
  g <- g + scale_color_discrete(name="Phenotype_CLOZUK.BGE")
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
  
  gglist[[i]] <- g
}

multiplot(plotlist = gglist, cols = 3)
screeplot(testing_experiments, main="Scree Plot", xlab="Components",ylim = c(0,1.4))
Testing_density_PC1 <- data.frame(data = testing_experiments$x[,1], as.factor(PCA_matrix_df$PHENOTYPE))
names(Testing_density_PC1) <- c("PC1","Phenotype")
ggplot(Testing_density_PC1, aes(x=PC1)) + geom_density(aes(group=Phenotype, colour = Phenotype, fill = Phenotype), alpha = 0.3)


for ( i in 1:15){
  
  g <- ggbiplot(testing2, obs.scale = 1, var.scale = 1, ellipse = F, choices = c(plots_index[1,i],plots_index[2,i]),
                circle = F , alpha = 1,varname.size = 3,groups = as.factor(PCA_matrix_df_cases$Colours))
  g <- g + scale_color_discrete(name="Phenotype_CLOZUK.BGE")
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
  
  gglist[[i]] <- g
}
multiplot(plotlist = gglist, cols = 3)

# PC1 various density plots
Testing_density_PC1 <- data.frame(data = testing$x[,1], as.factor(PCA_matrix_df$PHENOTYPE), as.factor(PCA_matrix_df$Colours))
names(Testing_density_PC1) <- c("PC1","Phenotype","Sample")

ggplot(Testing_density_PC1, aes(x=PC1)) + geom_density(aes(group=Phenotype, colour = Phenotype, fill = Phenotype), alpha = 0.3)
ggplot(Testing_density_PC1, aes(x=PC1)) + geom_density(aes(group=Sample, colour = Sample, fill = Sample), alpha = 0.3)

Testing_density_PC1_cases <- data.frame(data = testing2$x[,1], as.factor(PCA_matrix_df_cases$PHENOTYPE), as.factor(PCA_matrix_df_cases$Colours))
names(Testing_density_PC1_cases) <- c("PC1","Phenotype","Sample")

ggplot(Testing_density_PC1_cases, aes(x=PC1)) + geom_density(aes(group=Phenotype, colour = Phenotype, fill = Phenotype), alpha = 0.3)
ggplot(Testing_density_PC1_cases, aes(x=PC1)) + geom_density(aes(group=Sample, colour = Sample, fill = Sample), alpha = 0.3)



# PC2 Various density plots
Testing_density_PC2 <- data.frame(data = testing$x[,2], as.factor(PCA_matrix_df$PHENOTYPE), as.factor(PCA_matrix_df$Colours))
names(Testing_density_PC2) <- c("PC2","Phenotype", "Sample")

ggplot(Testing_density_PC2, aes(x=PC2)) + geom_density(aes(group=Phenotype, colour = Phenotype, fill = Phenotype), alpha = 0.3)
ggplot(Testing_density_PC2, aes(x=PC2)) + geom_density(aes(group=Sample, colour = Sample, fill = Sample), alpha = 0.3)

# PC2 cases density plot
Testing_density_PC2_cases <- data.frame(data = testing2$x[,2], as.factor(PCA_matrix_df_cases$PHENOTYPE), as.factor(PCA_matrix_df_cases$Colours))
names(Testing_density_PC2_cases) <- c("PC2","Phenotype", "Sample")

ggplot(Testing_density_PC2_cases, aes(x=PC2)) + geom_density(aes(group=Phenotype, colour = Phenotype, fill = Phenotype), alpha = 0.3)
ggplot(Testing_density_PC2_cases, aes(x=PC2)) + geom_density(aes(group=Sample, colour = Sample, fill = Sample), alpha = 0.3)


kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
p <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = plotly_POST(p, filename="surface/2")
chart_link

Testing_density_PC3 <- data.frame(data = testing$x[,3], as.factor(PCA_matrix_df$PHENOTYPE))
names(Testing_density_PC3) <- c("PC3","Phenotype")

ggplot(Testing_density_PC3, aes(x=PC3)) + geom_density(aes(group=Phenotype, colour = Phenotype, fill = Phenotype), alpha = 0.3)

Testing_density_PC4 <- data.frame(data = testing$x[,4], as.factor(PCA_matrix_df$PHENOTYPE))
names(Testing_density_PC4) <- c("PC4","Phenotype")

ggplot(Testing_density_PC4, aes(x=PC4)) + geom_density(aes(group=Phenotype, colour = Phenotype, fill = Phenotype), alpha = 0.3)

Testing_density_PC5 <- data.frame(data = testing$x[,5], as.factor(PCA_matrix_df$PHENOTYPE))
names(Testing_density_PC5) <- c("PC5","Phenotype")

ggplot(Testing_density_PC5, aes(x=PC5)) + geom_density(aes(group=Phenotype, colour = Phenotype, fill = Phenotype), alpha = 0.3)

Testing_density_PC6 <- data.frame(data = testing$x[,6], as.factor(PCA_matrix_df$PHENOTYPE))
names(Testing_density_PC6) <- c("PC6","Phenotype")

ggplot(Testing_density_PC6, aes(x=PC6)) + geom_density(aes(group=Phenotype, colour = Phenotype, fill = Phenotype), alpha = 0.3)


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
plot(pca_iris_rotated$scores[,4], pca_iris_rotated$scores[,1]), col = new_table$colour)
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
      col=summary(dat.em, data=dat)$classification, cex=2,pch='*')


library(plotly)
library(MASS)
# 3D Kernel Density distribution #

kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
p <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()


kd_PCA <- kde2d(Testing_density_PC1$PC1,Testing_density_PC2$PC2, n = 50)
kd_PCA_cases <- kde2d(Testing_density_PC1$PC1[Testing_density_PC1$Phenotype == 1],Testing_density_PC2$PC2[Testing_density_PC2$Phenotype==1], n = 50)
kd_PCA_controls <- kde2d(Testing_density_PC1$PC1[Testing_density_PC1$Phenotype == 0],Testing_density_PC2$PC2[Testing_density_PC2$Phenotype==0], n = 50)

kd_PCA_only_cases <- kde2d(Testing_density_PC1_cases$PC1, Testing_density_PC2_cases$PC2, n = 50)
cases_PCA <-plot_ly(x = kd_PCA_only_cases$x, y = kd_PCA_only_cases$y, z = kd_PCA_only_cases$z) %>% add_surface()

  combined_plot <- plot_ly() %>%
  add_surface(x = kd_PCA_cases$x, y = kd_PCA_cases$y, z = kd_PCA_cases$z, opacity = 1) %>%
  add_surface(x = kd_PCA_controls$x, y = kd_PCA_controls$y, z = kd_PCA_controls$z, opacity = 1)
  
# Cases plot
  theta <- seq(0,2*pi,length.out = 100)
  circle <- data.frame(x = cos(theta), y = sin(theta))
  p <- ggplot(circle,aes(x,y)) + geom_path()
  
  loadings <- data.frame(testing2$rotation, 
                         .names = row.names(testing2$rotation))
  p + geom_text(data=loadings, 
                mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) +
    coord_fixed(ratio=1) +
    labs(x = "PC1", y = "PC2")

# Controls plot
  theta <- seq(0,2*pi,length.out = 100)
  circle <- data.frame(x = cos(theta), y = sin(theta))
  p <- ggplot(circle,aes(x,y)) + geom_path()
  
  loadings <- data.frame(testing$rotation, 
                         .names = row.names(testing$rotation))
  p + geom_text(data=loadings, 
                mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) +
    coord_fixed(ratio=1) +
    labs(x = "PC1", y = "PC2")
  
  
p_PCA <- plot_ly(x = kd_PCA$x, y = kd_PCA$y, z = kd_PCA$z) %>% add_surface(surfacecolor = Testing_density_PC1$Phenotype)
p_PCA_cases <- plot_ly(x = kd_PCA_cases$x, y = kd_PCA_cases$y, z = kd_PCA_cases$z) %>% add_surface(surfacecolor = "red", opacity = 0.5, autocolorscale = F) %>%
plot_ly(x = kd_PCA_controls$x, y = kd_PCA_controls$y, z = kd_PCA_controls$z) %>% add_surface(surfacecolor = "blue", opacity = 0.5, auto)


p_PCA
# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = api_create(p, filename="surface/2")
chart_link

help(signup, package = 'plotly')



# Density plots of PRS cases and controls Schiz

Testing_density_PC1 <- data.frame(data = PCA_matrix_df$Schizophrenia, as.factor(PCA_matrix_df$PHENOTYPE), as.factor(PCA_matrix_df$Colours))
names(Testing_density_PC1) <- c("SCZ","Phenotype","Sample")

ggplot(Testing_density_PC1, aes(x=SCZ)) + geom_density(aes(group=Phenotype, colour = Phenotype, fill = Phenotype), alpha = 0.3)
