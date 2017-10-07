# ROC curves
library(RCurl)
x <- getURL("https://raw.githubusercontent.com/jbryer/CompStats/master/Data/titanic3.csv")
df <- read.csv(text = x)
log_reg <- function(df, size=10) {
  N <- nrow(df)
  size=10
  
  df <- df[sample(N),]
  
  num <- floor(N/size)
  rest <- N - num * size
  ncv <- cumsum(c(rep(size,num), rest))
  
  predictions <- data.frame(survived = df$survived, pred = NA)
  
  for(n in ncv) {
    v <- rep(TRUE, N)
    v[(n-size+1):n] <- FALSE
    
    lr <- glm(survived ~ ., data = df[v,], family = binomial(logit))
    predictions[!v,"pred"] <- predict(lr, newdata=df[!v,], type="response")
  }
  
  return(predictions)
}

predictions <- log_reg(df, size = 10)

str(df)


category <- c(1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0)
prediction <- rev(seq_along(category))
prediction[9:10] <- mean(prediction[9:10])

library(pROC)
roc_obj <- roc(category, prediction)
auc(roc_obj)

roc_df <- data.frame(
  TPR=rev(roc_obj$sensitivities), 
  FPR=rev(1 - roc_obj$specificities), 
  labels=roc_obj$response, 
  scores=roc_obj$predictor)

testing <- prcomp(PCA_matrix_df[c(3:7)],center = T)


roc_obj1 <- roc(PCA_matrix_df$PHENOTYPE, testing$x[,1], plot = T, ci = T)
roc_obj2 <- roc(PCA_matrix_df$PHENOTYPE, testing$x[,2], plot = T, ci = T)
roc_obj3 <- roc(PCA_matrix_df$PHENOTYPE, testing$x[,3], plot = T, ci = T)
roc_obj4 <- roc(PCA_matrix_df$PHENOTYPE, testing$x[,4], plot = T, ci = T)
roc_obj5 <- roc(PCA_matrix_df$PHENOTYPE, testing$x[,5], plot = T, ci = T)
roc_obj6 <- roc(PCA_matrix_df$PHENOTYPE, testing$x[,6], plot = T, ci = T)


roc_obj1 <- roc(PCA_matrix_2$PHENOTYPE, testing_reduce_controls$x[,1], plot = T, ci = T)
roc_obj_limited_controls <- roc(PCA_matrix_2$PHENOTYPE, testing_reduce_controls$x[,1], plot = T, ci = T)
roc_obj2 <- roc(PCA_matrix_2$PHENOTYPE, testing_reduce_controls$x[,2], plot = T, ci = T)
roc_obj3 <- roc(PCA_matrix_2$PHENOTYPE, testing_reduce_controls$x[,3], plot = T, ci = T)
roc_obj4 <- roc(PCA_matrix_2$PHENOTYPE, testing_reduce_controls$x[,4], plot = T, ci = T)
roc_obj5 <- roc(PCA_matrix_2$PHENOTYPE, testing_reduce_controls$x[,5], plot = T, ci = T)
roc_obj6 <- roc(PCA_matrix_2$PHENOTYPE, testing_reduce_controls$x[,6], plot = T, ci = T)

roc_obj_SCZ <- roc(PCA_matrix_2$PHENOTYPE, PCA_matrix_2$Schizophrenia, plot = T, ci = T)

roc_obj_SCZ <- roc(PCA_matrix_2$PHENOTYPE, PCA_matrix_2$Schizophrenia, plot = T, ci = T)
roc_obj_BIP <- roc(PCA_matrix_2$PHENOTYPE, PCA_matrix_2$Bipolar, plot = T, ci = T)
roc_obj_EDU <- roc(PCA_matrix_2$PHENOTYPE, PCA_matrix_2$Educational_attainment, plot = T, ci = T)
roc_obj_NEU <- roc(PCA_matrix_2$PHENOTYPE, PCA_matrix_2$Neuroticism, plot = T, ci = T)
roc_obj_BIPvsSCZ <- roc(PCA_matrix_2$PHENOTYPE, PCA_matrix_2$BIPvsSCZ, plot = T, ci = T)
roc_obj_MDD <- roc(PCA_matrix_2$PHENOTYPE, PCA_matrix_2$PGC_MDD, plot = T, ci = T)
roc_obj_IQ <- roc(PCA_matrix_2$PHENOTYPE, PCA_matrix_2$IQ, plot = T, ci = T)


roc_obj_wo_SCZ <- roc(PCA_matrix_2$PHENOTYPE, testing_wo_SCZ$x[,1], plot = T, ci = T)
roc_obj_wo_BIP <- roc(PCA_matrix_2$PHENOTYPE, testing_wo_BIP$x[,1], plot = T, ci = T)
roc_obj_wo_EDU <- roc(PCA_matrix_2$PHENOTYPE, testing_wo_EDU$x[,1], plot = T, ci = T)
roc_obj_wo_MDD <- roc(PCA_matrix_2$PHENOTYPE, testing_wo_MDD$x[,1], plot = T, ci = T)
roc_obj_wo_BIPvsSCZ <- roc(PCA_matrix_2$PHENOTYPE, testing_wo_BIPvsSCZ$x[,1], plot = T, ci = T)
roc_obj_wo_Neurot <- roc(PCA_matrix_2$PHENOTYPE, testing_wo_Neurot$x[,1], plot = T, ci = T)
roc_obj_wo_IQ <- roc(PCA_matrix_2$PHENOTYPE, testing_wo_IQ$x[,1], plot = T, ci = T)
roc_obj_wo_Intell <- roc(PCA_matrix_2$PHENOTYPE, testing_wo_intelligence$x[,1], plot = T, ci = T)

auc(roc_obj_wo_SCZ)
roc_obj_wo_SCZ$ci
auc(roc_obj_wo_BIP)
roc_obj_wo_BIP$ci
auc(roc_obj_wo_EDU)
roc_obj_wo_EDU$ci
auc(roc_obj_wo_MDD)
roc_obj_wo_MDD$ci
auc(roc_obj_wo_BIPvsSCZ)
roc_obj_wo_BIPvsSCZ$ci
auc(roc_obj_wo_Neurot)
roc_obj_wo_Neurot$ci



roc_obj_exper <- roc(PCA_matrix_df$PHENOTYPE, testing_experiments$x[,1])
roc_obj_test <- roc(PCA_matrix_df$PHENOTYPE, PCA_matrix_df$Schizophrenia, plot=T)

auc(roc_obj1)
roc_obj1$ci
auc(roc_obj2)
roc_obj2$ci
auc(roc_obj3)
roc_obj3$ci
auc(roc_obj4)
roc_obj4$ci
auc(roc_obj5)
roc_obj5$ci
auc(roc_obj6)
roc_obj6$ci


auc(roc_obj_SCZ)
roc_obj_SCZ$ci
auc(roc_obj_BIP)
roc_obj_BIP$ci
auc(roc_obj_EDU)
roc_obj_EDU$ci
auc(roc_obj_NEU)
roc_obj_NEU$ci
auc(roc_obj_BIPvsSCZ)
roc_obj_BIPvsSCZ$ci
auc(roc_obj_MDD)
roc_obj_MDD$ci


