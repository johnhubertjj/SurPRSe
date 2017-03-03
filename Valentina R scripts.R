#### Valentina scripts ####

# Extracting common SNPs #

Data:  /home/SHARED/forJudith


Extract common SNPs

dat<-read.table(file="xxxxxx.summary", header=T)
o<-order(dat$BP)
dat<-dat[o,]
o<-order(dat$CHR)
dat<-dat[o,]
dat$LOC<-paste(dat$CHR,":",dat$BP, sep="")

thy<-read.table(file="xxxxx.bim", header=F)
thy$LOC<-paste(thy$V1,":",thy$V4, sep="")

m<-merge(dat, thy, by.x="LOC", by.y="LOC",  sort=F, all=F)
cat("common N=",nrow(m),"\n")




save the list: common SNPs: CHE:POS

clump r2=0.2, p1 0.5, p2 0.5

# Merging 2 SNP tables and flipping SNP's #

m<-merge(clump, data, by.x="LOC", by.y="LOC", all=F, sort=F)

m$A1<-toupper(m$A1)
m$A2<-toupper(m$A2)

cat(fout, nrow(m))
a<-which(m$A1=="C" & m$A2=="G"); if (length(a)>0) {m<-m[-a,]}
a<-which(m$A1=="G" & m$A2=="C"); if (length(a)>0) {m<-m[-a,]}
a<-which(m$A1=="A" & m$A2=="T"); if (length(a)>0) {m<-m[-a,]}
a<-which(m$A1=="T" & m$A2=="A"); if (length(a)>0) {m<-m[-a,]}

cat(fout, ", remove A-T, C-G: N=", nrow(m))


#FLIPPING ---------------------------------------------------
a<-which(m$V5==m$A2 & m$V6==m$A1)
m$B[a]<- (-m$B[a])
m$A1[a]<-as.character(m$V5[a])
m$A2[a]<-as.character(m$V6[a])



a<-which(m$V5==m$A1 & m$V6==m$A2)
d<-seq(1:nrow(m));d<-d[-a]

b1<-which(m$A1[d]=="A")
b2<-which(m$A1[d]=="C")
b3<-which(m$A1[d]=="G")
b4<-which(m$A1[d]=="T")
m$A1[d[b1]]<-"T"
m$A1[d[b2]]<-"G"
m$A1[d[b3]]<-"C"
m$A1[d[b4]]<-"A"

b1<-which(m$A2[d]=="A")
b2<-which(m$A2[d]=="C")
b3<-which(m$A2[d]=="G")
b4<-which(m$A2[d]=="T")
m$A2[d[b1]]<-"T"
m$A2[d[b2]]<-"G"
m$A2[d[b3]]<-"C"
m$A2[d[b4]]<-"A"

a<-which(m$V5==m$A2 & m$V6==m$A1)
m$B[a]<- (-m$B[a])
m$A1[a]<-as.character(m$V5[a])
m$A2[a]<-as.character(m$V6[a])

a<-which(m$V5!=m$A1 | m$V6!=m$A2)
if (length(a)>0) m<-m[-a,]

cat(fout, ", remove A-T, C-G and other mismatches: N=", nrow(m))

#### SCORING:------------------------------------------------------------------------
  
  sig<-c(1e-4,1e-3,1e-2,0.05,0.1,0.2,0.3,0.4,0.5)

fout<-"XXXXXX.txt"
fcl<-"clump_r0.2_1000kb.clumped.snps"
fname<-"PGC_noCLOZUKSummary.txt"

clump<-read.table(file=fcl, header=F)
data<-read.table(file=fname, header=T)
m<-merge(clump, data, by.x="V1", by.y="SNP", all=F, sort=F)
m$Allele1<-toupper(m$Allele1)
m$Allele2<-toupper(m$Allele2)

#score
for (i in 1:(length(sig)))
{
  fn<-sprintf("scores/%s_%f.score", fout, sig[i])
  a<-which(m$P.value<=sig[i])
  if (length(a)>0)
  {
    tmp<-m[a,]
    out<-tmp[,c("V1","Allele1","Effect")]
  }
  if (nrow(out)>0) {write.table(file=fn, out, row.names=F, col.names=F, quote=F, sep="\t")}
}





#### PROFILING: ------------------------------------------------------------------------

sig<-c(1e-4,1e-3,1e-2,0.05,0.1,0.2,0.3,0.4,0.5)

for (i in 1:(length(sig)))
{
  for (chr in 1:22)
  {
    command_txt<-sprintf("plink2 --silent --bfile  XXX_chr%i --score scores/%s_%f.score,  --out  profiles/chr%i.clump_r0.2_1000kb_%f", chr, sig[i], chr, sig[i])
    
    system(command_txt)
  }
} #sig




##### LOGISTIC REGRESSION FOR EACH GENE:  -----------------------------------------------
  
  Cov<-CLOZUK2.r5_Final.eigenvec

for (i in 1:length(sig))
{
  profn<-sprintf("clumped.maf0.1_rep_common.out_%f.profile", sig[i])
  dat<-read.table(file=profn, header=T)
  dat<-merge(cov,dat,by.x="FID", by.y="FID", all=F)
  
  res$model[i]<-sig[i]
  
  model0<-glm(SCORE~PC1+PC2+PC3+PC4+PC5+PC6+PC8+PC11+PC12+PC13+PC16+PC18,family=gaussian, data=dat)
  m1<-mean(residuals(model0))
  sd1<-sd(residuals(model0))
  dat$NORMSCORE<-(residuals(model0)-m1)/sd1
  #             hist(dat$NORMSCORE)
  model<-glm(PHENO~dat$NORMSCORE, family=binomial, data=dat)
  
  #             res$Effect[i]<-summary(model)$coefficients[2, "Estimate"]
  #             res$SE[i]<-summary(model)$coefficients[2, "Std. Error"]
  res$p[i]<- summary(model)$coefficients[2, "Pr(>|z|)"]
  #             res$NSNPs[i]<-max(dat$CNT)/2
}#i
mtext("FT4-Replication. MAF>0.1", side = 1, line = -1, outer = TRUE)
write.table(file=fout, res, col.names=T, row.names=F, quote=F, sep="\t")
