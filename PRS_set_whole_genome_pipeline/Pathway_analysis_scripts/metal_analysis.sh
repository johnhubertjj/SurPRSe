# This file includes a series of comments. Each comment is marked 
# by a # character as the first character in the line.
#
# This is a comment!

# Meta-analysis weighted by standard error does not work well
# when different studies used very different transformations.
# In this case, some attempt was made to use similar trait
# transformation and you can request a standard error based
# analysis by uncommenting the following line:
# SCHEME   STDERR

# Usually, it is a good to apply genomic control to each
# input file. However, in this example, all the markers being
# examined reside in strongly associated loci and it probably 
# is not a good idea. To find out what happens if you turn
# on genomic control, uncomment the following line.
# GENOMICCONTROL ON

# To help identify allele flips, it can be useful to track
# allele frequencies in the meta-analysis. To enable this
# capability, uncomment the following two lines.
# AVERAGEFREQ ON
# MINMAXFREQ ON

# To restric meta-analysis to two previously reported SNPs
# and summarize study specific results, uncomment the two
# lines that follow.
# ADDFILTER SNP IN (rs10830963,rs563694)
# VERBOSE ON
# Describe and process the DGI input files
MARKER   RS_SNP
ALLELE   A1 A2
#DEFAULT 69416
EFFECT   log(OR)
SCHEME STDERR
STDERR   SE
PVAL     P
VERBOSE ON
ADDFILTER INFO > 0.9
PROCESS /home/johnhubert/Documents/CLOZUK_cognition_METAL/CLOZUK_GWAS-noCOGS/PGC2noCLOZUK_whole_genome_flipped_RS_SNP.txt
#DEFAULT 4173
PROCESS /home/johnhubert/Documents/CLOZUK_cognition_METAL/CLOZUK_GWAS-noCOGS/CLOZUK2_noPGC2_whole_genome_flipped_RS_SNP.txt

OUTFILE CLOZUK_PGC2_noCOGS .tbl
# Execute meta-analysis
ANALYZE

