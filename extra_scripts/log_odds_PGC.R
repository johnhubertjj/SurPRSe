## prerequisites
library(data.table)

##setting working directory
wd<-getwd()
#library(parallel)

## reading in file 
d = fread(paste(wd,"/daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt",sep = ""))

## converting odds ratio
odds_ratio_row<-d[,OR]
BETA<-log(odds_ratio_row)

## adding BETA column  
New_PGC<-cbind(beta_PGC_summary_stats,beta2)
New_PGC<- New_PGC[,c(1:8,18,9:17)]
write(BETA,file=paste(wd,"/PGC_summ_stats_beta_coeff.txt"))
write(New_PGC,file=paste(wd,"/New_PGC_result_file.txt"))

#The rest will work on a cluster, but not on desktop##
#Probably better as a python script

#system.time(log(odds_ratio_row))


#no_cores <- detectCores() - 1
#cl <- makeCluster(no_cores)

#test1<-parSapply(cl, odds_ratio_row,function(odds_ratio){ log(odds_ratio) })






#stopCluster(cl)
