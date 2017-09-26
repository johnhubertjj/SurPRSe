library(data.table)
require(cowplot)

PRS.profiles

Bipolar1 <-fread("~/Dropbox/PRS_scores_different_summary_stats_datasets/BIP_CLOZUK_whole_genome_significance_threshold_at_0.01.profile")
Bipolar1$SCORE <- scale(Bipolar1$SCORE,center = T,scale = T)
Bipolar2 <-fread("~/Dropbox/PRS_scores_different_summary_stats_datasets/BIP_CLOZUK_whole_genome_significance_threshold_at_0.25.profile")
Bipolar2$SCORE <- scale(Bipolar2$SCORE,center = T,scale = T)
Bipolar3 <-fread("~/Dropbox/PRS_scores_different_summary_stats_datasets/BIP_CLOZUK_whole_genome_significance_threshold_at_0.5.profile")
Bipolar3$SCORE <- scale(Bipolar3$SCORE,center = T,scale = T)

ggplot(data=Bipolar1, aes(Bipolar1$SCORE)) + 
  geom_histogram( 
                 col="red", 
                 aes(fill=..count..))

gc <- ggplot(data=Bipolar2, aes(Bipolar2$SCORE)) 
gc <- gc +  geom_histogram( 
            col="red", 
            aes(fill=..count..))
gc <- gc + theme(axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line("black"),
                 legend.title = element_blank(),
                 legend.position = 'none'
                 )
gc

ga <- ggplot(data=Bipolar1, aes(Bipolar1$SCORE)) 
ga <- ga +  geom_histogram( 
  col="red", 
  aes(fill=..count..))
ga <- ga + theme(axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line("black"),
                 legend.title = element_blank(),
                 legend.position = 'none'
)
ga

gs <- ggplot(data=Bipolar3, aes(Bipolar3$SCORE)) 
gs <- gs +  geom_histogram( 
  col="red", 
  aes(fill=..count..))
gs <- gs + theme(axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line("black"),
                 legend.title = element_blank(),
                 legend.position = 'none'
)


#plot_grid(gc,ga,gs,labels=c("0.01","0.25","0.5"),align = 'h',nrow = 1)
plot_grid(gc,ga,gs,align = 'h',nrow = 1)
