## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
library(GSA)

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("mygene")
library(mygene)

Calcium_ion_ID <- read.table("~/Dropbox/GO_0070509.txt",header = F, stringsAsFactors = F,fill = T)
Synaptic_transmission <- read.table("~/Dropbox/GO_0007268.txt", header=F,stringsAsFactors = F,fill = T)
Membrane_depolarization_during_action_potential <- read.table("~/Dropbox/GO_0086010.txt", header=F,stringsAsFactors = F,fill = T)
  
Calcium_test1 <- Calcium_ion_ID$V2
Calcium_test1 <- unique(Calcium_test1)

Membrane_depolr_test1 <- Membrane_depolarization_during_action_potential$V2
Membrane_depolr_test1 <- unique(Membrane_depolr_test1)

Synaptic_transmission_test1 <- Synaptic_transmission$V2
Synaptic_transmission_test1 <- unique(Synaptic_transmission_test1)


a <- queryMany(Calcium_test1, scopes="symbol", fields="entrezgene", species="human",)
b <- queryMany(Membrane_depolr_test1, scopes="symbol", fields="entrezgene", species="human")
c <- queryMany(Synaptic_transmission_test1, scopes="symbol", fields="entrezgene", species="human",returnall=T)

a$Pathway_name <- rep("Calcium_ion_import_GO0070509",nrow(a))

b$Pathway_name <- rep("Membrane_depolarization_during_action_potential_GO0086010",nrow(b))
c$response$Pathway_name <- rep("Synaptic_transmission_GO0007268",nrow(c$response))

Pocklington <- fread("~/Dropbox/Stationary_data/Pocklington2015_134sets_LoFi.txt")
setkey(Pocklington,V1)
#Neuropsychiatric_datasets[[i]] <- Neuropsychiatric_datasets[[i]][grep(paste(Groups_to_keep,collapse="|"), 
#Neuropsychiatric_datasets[[i]]$FID, value=TRUE)]

Pathways_from_pocklington <- c("\\babnormal_behavior\\b", "\\b5HT_2C\\b", "\\bFMRP_targets\\b", "\\babnormal_nervous_system_electrophysiology\\b", "\\bCav2_channels\\b", "\\babnormal_long_term_potentiation\\b", "\\babnormal_grooming_behavior\\b", "\\bLek2015_LoFintolerant_90\\b")

Pathways_to_use <- Pocklington[grepl(paste(Pathways_from_pocklington, collapse = "|"),Pocklington$V1)]
Other_pathways <- rbind(a[,c("Pathway_name","entrezgene")], b[,c("Pathway_name","entrezgene")], c$response[,c("Pathway_name","entrezgene")])
Other_pathways2 <- Other_pathways[!is.na(Other_pathways$entrezgene),]
colnames(Other_pathways2) <- c("V1","V2")
All_pathways <- rbind(Pathways_to_use, Other_pathways2)

write.table(Pathways_to_use, file = "~/Dropbox/Selected_Pocklington_pathways.txt",row.names = F, quote = F, col.names = F)
write.table(Other_pathways2, file = "~/Dropbox/Antonio_GO_enriched_pathways.txt",row.names = F, quote = F, col.names = F)
write.table(All_pathways, file = "~/Dropbox/Selected_Pocklington_plus_GO_pathways_SCHIZ.txt",row.names = F, quote = F, col.names = F)


write.table(Pathw)
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
a <-getGene(test1, type="entrezgene", mart=mart)

enrichr_GO <- GSA.read.gmt("~/Dropbox/GO_Biological_Process_2017.txt")
