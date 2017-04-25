###################
# MANDATORY STEPS #
###################

# Set Working directory (project should start in dbSNP directory)
setwd("./dbSNP")

# Enable JIT for better performance on hidden loops
library(compiler)
enableJIT(3)

#############################
# Load and update BIM files #
#############################

# Load dbSNPb142 (GRCh37 version)
library("data.table")
load("/Users/johnhubert/Dropbox/RS_dbSNP_script/dbSNP/dbSNP_GRCh37p13_db149.RData")

# Function to update the RS IDs (now with 100% more progress bars!)
library(tcltk2)
updateRS <- function(bimframe) {
  pb <- tkProgressBar(title = "Progress Bar for UpdateRS", min = 0,max = nrow(bimframe), width = 300)
  newbim <- bimframe 
  rsvec <- as.vector(bimframe$RS)
  chrvec <- as.vector(bimframe$CHR)
  bpvec <- as.vector(bimframe$BP) 
  Avec <- as.vector(bimframe$A)
  Bvec <- as.vector(bimframe$B)
  for (i in 1:nrow(bimframe)) {
    if (chrvec[i] == 0) {rsvec[i] <- paste0(chrvec[i],":",rsvec[i])} # special rename for CHR0 variants (will be discarded anyway)
    else {
      #if (length(grep("rs",rsvec[i])) == 1) {next} # Skip variant if this is already named
      rsid <- eval(parse(text=paste0("chr",chrvec[i])))[list(bpvec[i])]$RS[1]
      if (is.na(rsid)) {rsvec[i] <- paste0(chrvec[i],":",bpvec[i])} # If variant is not found, use 1KGP naming convention minus the REF/ALT alleles (as SNPs might not be in plus strand)
      else {rsvec[i] <- paste0("rs",rsid) }
      setTkProgressBar(pb, i,label=paste(i,"in",nrow(bimframe),"SNPs processed"))}}
  newbim$RS <- rsvec
  close(pb)
  return(newbim)
}

# Function to select duplicated SNPs at random (for name changing)
duplicated.random = function(x, incomparables = FALSE, ...)
{
  if ( is.vector(x) )
  {
    permutation = sample(length(x))
    x.perm      = x[permutation]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
  else if ( is.matrix(x) )
  {
    permutation = sample(nrow(x))
    x.perm      = x[permutation,]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
  else
  {
    stop(paste("duplicated.random() only supports vectors",
               "matrices for now."))
  }
}


############################
# PROCESS CLOZUK2 DATASETS #
############################

# Generation Scotland data (use a BIM.BAK file to preserve original)
#GenSCOT <- fread("GenScot-qc.bim.bak")
GenSCOT <- fread("BPSCZ.bp_v_scz.results_GR37.p13.txt")

setnames(GenSCOT, c("CHR","RS","BP","A","B","INFO","OR","SE","P", "ngt"))


GenSCOT <- updateRS(GenSCOT)

# Remove duplicates (if any)
GenSCOT.rnddup <- which(duplicated.random(GenSCOT$RS))
GenSCOT[GenSCOT.rnddup,]$RS <- paste0(GenSCOT[GenSCOT.rnddup,]$RS,".dup")
GenSCOT.duplicates <- as.data.frame(GenSCOT[grep(".dup",GenSCOT$RS),]$RS)
GenSCOT.rnddup <- which(duplicated.random(GenSCOT$RS))
GenSCOT[GenSCOT.rnddup,]$RS <- paste0(GenSCOT[GenSCOT.rnddup,]$RS,".2")
GenSCOT.duplicates <- as.data.frame(GenSCOT[grep(".dup",GenSCOT$RS),]$RS)

# Write new result
write.table(GenSCOT,file="using_rs_iD_converter_hg19_BIZvsSCZ.txt",quote = F,row.names = F,col.names=F)
write.table(GenSCOT.duplicates,file="GenScot-qc.dup",quote = F,row.names = F,col.names=F)
rm(GenSCOT)
rm(GenSCOT.duplicates)
rm(GenSCOT.rnddup)

map_file_for_liftover <- GenSCOT[,.(CHR, SNP, 0, BP)]
Bed_file_for_liftover <- GenSCOT[,BP0 := (BP-1)]

Bed_file_for_liftover$CHR <- sub("^", "chr", Bed_file_for_liftover$CHR)
Bed_file_for_liftover <- Bed_file_for_liftover[,.(CHR,BP0, BP, SNP, A1, A2)]
setcolorder(Bed_file_for_liftover, c(1,11,3,2,4:10))
write.table(Bed_file_for_liftover, file = "BIPvsSCZ_bed_file_for_liftover.txt", quote = F, row.names = F, col.names = F)
write.table(map_file_for_liftover, file = "BIPvsSCZ_map_file_for_liftover.txt", quote = F, row.names = F, col.names = F)


new_bp_positions <- fread("~/Downloads/remapped_BIPvsSCZ_map_file_for_liftover.txt.bed")
setnames(new_bp_positions, c("CHR", "BP_dont_use", "BP_use", "SNP"))
new_bp_positions$CHR <- sub("chr", "", new_bp_positions$CHR)
new_bp_positions$CHR <- as.integer(new_bp_positions$CHR)
a <- merge(new_bp_positions, GenSCOT, by = c("SNP", "CHR") , all = F)
