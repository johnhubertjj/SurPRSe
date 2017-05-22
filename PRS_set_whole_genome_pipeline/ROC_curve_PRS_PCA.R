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

roc_obj2 <- roc(PCA_matrix_df$PHENOTYPE, testing$x[,1])
roc_obj3 <- roc(PCA_matrix_df$PHENOTYPE, testing_experiments$x[,1])
roc_obj4 <- roc(PCA_matrix_df$PHENOTYPE, PCA_matrix_df$Schizophrenia)
auc(roc_obj2)
auc(roc_obj3)
auc(roc_obj4)
