library(data.table)
Pocklington_pathways <- read.table("~/Dropbox/Stationary_data/Pocklington2015_134sets_LoFi.txt",stringsAsFactors = F)
Pocklington_morphology_2 <- Pocklington_pathways[grepl("morphology",Pocklington_pathways$V1) | 
                                                   grepl("thin_cerebral_cortex",Pocklington_pathways$V1) |
                                                   grepl("abnormal_brain_size",Pocklington_pathways$V1),]

Pocklington_morphology_2$V1 <- "Morphology_schizophrenia_related"

Pocklington_morphology_2_unique<- Pocklington_morphology_2[!(duplicated(Pocklington_morphology_2$V2)),]

Pocklington_without_morphology <- Pocklington_pathways[!(grepl("morphology",Pocklington_pathways$V1) | 
                                                     grepl("thin_cerebral_cortex",Pocklington_pathways$V1) |
                                                     grepl("abnormal_brain_size",Pocklington_pathways$V1)),]

Pocklington_without_morphology$V1 <- "MGI_Schizophrenia_related_excluding_morphology_including_LoF"

Pocklington_without_morphology_deduplicated<- Pocklington_without_morphology[!(duplicated(Pocklington_without_morphology$V2)),]

Pocklington_changed_terms_names <- rbind(Pocklington_without_morphology_deduplicated,Pocklington_morphology_2_unique)

write.table(Pocklington_morphology_2_unique,file = "~/Dropbox/Stationary_data/Pocklington2015_134sets_LoFi_morphology_only_deduplicated.txt",row.names = F,col.names = F,quote = F)
write.table(Pocklington_changed_terms_names,file = "~/Dropbox/Stationary_data/Pocklington2015_134sets_LoFi_2sets_morphology_notmorphology_deduplicated.txt",row.names = F,col.names = F,quote = F)
