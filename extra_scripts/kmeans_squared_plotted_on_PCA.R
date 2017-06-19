set.seed(20)
kmeans_pc12 <- as.data.frame(testing2$x[,1:2])
kmeans_all_pcs <-as.data.frame(testing2$x)
testing_kmeans <- kmeans(kmeans_pc12, 50, nstart = 20)
testing_kmeans
testing_kmeans$cluster <- as.factor(testing_kmeans$cluster)

ggplot(kmeans_pc12, aes(PC1, PC2, color = testing_kmeans$cluster)) + geom_point()

kmeans_original <- as.data.frame(PCA_matrix_df_cases[,2:6])
testing_kmeans_original <- kmeans(kmeans_original, 6, nstart = 20)
testing_kmeans_original$cluster <- as.factor(testing_kmeans_original$cluster)
ggplot(kmeans_pc12, aes(PC1,PC2 , color = testing_kmeans_original$cluster)) + geom_point()
ggplot(kmeans_all_pcs, aes(PC1,PC3 , color = testing_kmeans_original$cluster)) + geom_point()
