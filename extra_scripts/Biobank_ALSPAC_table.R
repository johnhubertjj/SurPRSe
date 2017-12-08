library(ggplot2)

Training_data <- data.frame(
  Dataset = factor(c("PGC1","PGC1","PGC1swe","PGC1swe","PGC2","PGC2","CLOZUK","CLOZUK"), levels=c("PGC1","PGC1swe","PGC2","CLOZUK")),
  cohort = factor(c("Cases", "Controls","Cases", "Controls","Cases", "Controls","Cases", "Controls"), levels = c("Cases", "Controls")),
  Number_of_individuals = c(8831, 12067, 13833, 18310, 29415, 40101, 40675,64643)
)

p <- ggplot(data=Training_data, aes(x=Dataset, y=Number_of_individuals, fill=cohort)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_fill_manual(values=c("#999999", "#E69F00"))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
p <- p+ ggtitle("Training datasets")
p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))