
data(iris)
library(data.table)
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

PCA_matrix <- data.frame(data = merge1$FID,merge1$PHENO.y, merge1$NORMSCORE, merge2$NORMSCORE, merge3$NORMSCORE, merge4$NORMSCORE, merge5$NORMSCORE, merge6$NORMSCORE)
library(FactoMineR)
names(PCA_matrix) <- c("Individuals","PHENOTYPE", "Schizophrenia", "Bipolar", "Educational_attainment", "PGC_MDD", "BIPcsSCZ", "Neuroticism")

write.csv(PCA_matrix, file = "E:/PRS_cross_disorder_table.csv", col.names = T, row.names = F)
testing <- prcomp(PCA_matrix[,3:8],center = T)
biplot(testing)

a <- cov(PCA_matrix[2:8])
eigenPRS <- eigen(a)

plot(testing$ind$coord[,2], testing$ind$coord[,3])
a <- merge1$NORMSCORE 
b <- merge2$NORMSCORE
c <- merge3$NORMSCORE
d <- merge4$NORMSCORE
e <- merge5$NORMSCORE
f <- merge6$NORMSCORE

new.data.frame <- data.frame(PRS.Profiles.with.covariates$FID,a,b,c,d,e)

PCA <- prcomp(new.data.frame)
write.csv <- n

new_data_frame <- 
test_PCA <- merge(merge1,merge3,by = c("FID", "PHENO.y"), all = T)
new_table <- test_PCA[,.(PHENO.y,NORMSCORE.x,NORMSCORE.y)]
new_table2 <- test_PCA[PHENO.y==1,.(NORMSCORE.x,NORMSCORE.y)]
pca2 <- prcomp(new_table2,center = T,scale. = T)
plot(pca2)
Phenotype <- as.factor(new_table$PHENO.y)

library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(testing_wah, obs.scale = 1, var.scale = 1, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)


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

rawLoadings_cases <- my.prc_cases$rotation[,1:6] %*% diag(my.prc_cases$sdev, 6, 6)
rotatedLoadings <- varimax(rawLoadings_cases)$loadings
invloadings <- t(pracma::pinv(rotatedLoadings))
scores <- scale(PCA_matrix_cases[,-1:-2]) %*% invloadings


iris

my.abs     <- abs(cor(PCA_matrix_cases[,-1:-2]))
my.colors  <- dmat.color(my.abs)
my.ordered <- order.single(cor(PCA_matrix_cases[,-1:-2]))
cpairs(PCA_matrix_cases[,-1:-2], my.ordered, panel.colors=my.colors, gap=0.5)

my.abs     <- abs(cor(PCA_matrix_controls[,-1:-2]))
my.colors  <- dmat.color(my.abs)
my.ordered <- order.single(cor(PCA_matrix_controls[,-1:-2]))
cpairs(PCA_matrix_controls[,-1:-2], my.ordered, panel.colors=my.colors, gap=0.5)
