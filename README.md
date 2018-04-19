# Super-PRS

* A workflow in Bash, R and python designed to allow easy processing and production of polygenic risk scores on [ARCCA's](http://www.cardiff.ac.uk/advanced-research-computing) Raven supercomputer
* If used for the creation of results to be published (one can dream...), please acknowledge  John hubert <a href="mailto:hubertjj@cardiff.ac.uk">mailto:hubertjj@cardiff.ac.uk</a> and Valentina Escott-Price <a href="mailto:escottpricev@cardiff.ac.uk">mailto:escottpricev@cardiff.ac.uk</a>

# Current datasets tested in PRS analysis:

PGCnoCLOZUK
-------------- 
Location on rocks: /home/SHARED/PGC/daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt
Summary stats

CLOZUK_no_PGC:
---------------
Location on rocks: /mnt/databank/CLOZUK/GWAS/BGE/\*CLOZUK_GWAS_BGE\*.tar.gz
Build: 
Best guess genotype data


# Preparing the config file
While not ideal, all preparations for the pipeline should be completed within _PRS\_arguments.sh_ and _PRS\_Pipeline\_for\_local.sh_

in _PRS\_arguments.sh_ where the line:
> elif [ "$whereami" = 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then

is found, replace **v1711-0ab8c3db.mobile.cf.ac.uk** with the character string found after running following command in a terminal:
> uname -n

for the home_OS argument, add the physical location of your hard-drive where the the analysis will be run to replace the root "**/**".
For example in MAC, most users have '/Users' while Unix/windows usually have '/home'. if unsure, open a terminal session and type:
> pwd

and extract the string of the first directory after the root '**/**'

## Arguments

path\_to\_stationary\_data: Location of standardised files for input into your Polygenic risk score analysis: please use in the following format (for now:)

> ${home_OS}/path/to/stationary/data 

 path\_to\_covariate\_file: Location of population covariates to adjust Polygenic risk score if required. Text-delimited format with columns: FID IID PC1 PC2 etc... 

> ${path\_to\_stationary\_data}**covarate\_file**

**path\_to\_chromosome\_length\_file**: **FIND OUT WHAT THIS DOES IN THE PIPELINE**

**path\_to\_new\_fam\_file**: Location of alternative fam file in case of plink files from meta-analysis or use of an alternatve phenotype: please use the format above replacing the covariate file with plink fam file 


**training\_set\_name**: User defined name for the **Training** dataset in the polygenic risk score. Does not have to match the file name but **MUST NOT** include any delimiters.  



**path\_to\_training\_dataset**: path to the **Training dataset** for input into polygenic risk score.  


**validation\_set\_name**: User defined name for the **Test** dataset in the polygenic risk score. Does not have to match the file name but **MUST NOT** include any delimiters.  

**path\_to\_validation\_dataset**: path and plink prefix for **Test dataset**. File format must be the same as input into the --bfile argument to plink.

* The file format can be one set of plink files (aka one each of .bim/.bed/.fam)
 * **NOTE** all analysis is split by chromosome, please set the split\_by\_chromosome\_required argument to **TRUE** 
* Also accepts plink files separated into chromosomes in normal plink format or compressed (ONLY .tar.gz)
 * If compressed, the compressed volume and the prefix of the file names within **MUST** match this argument      
 * **NOTE** if this is the case, the split\_by\_chromosome\_required argument **MUST BE** set to **FALSE**
 * additionally, each file must have the character string "\_chr" followed by the chromosome posistion appended to validation\_set\_full\_name\_without\_chromosome
  * eg: dataset\_chr1.tar.gz

**Pathway\_filename**: path to gene-set file describing the locations of genes within gene-sets. Must be in format: gene-set gene

* For each gene, the specified gene-set name is repeated on each row until the gene-set is complete
 * Gene set names would be preferable without any special characters except underscores for delimiters

**Gene\_location\_filename**: path to Gene location file. Either in MAGMA 1.6 .gene.loc fileformat or with the following columns: ENTREZ\_GENE\_ID, CHR, BP\_START, BP\_END, STRAND, GENE\_NAME. Ensure extension to file is also included and the BPSTART and BPEND are of gene boundaries if you wish to utilise extended gene boundaries (see Gene\_regions argument)  
 
**split\_by\_chromosome\_required**: if plink files are already split by chromosome, please set to FALSE (without quotes), if one set of plink files is provided, please set to TRUE (also without quotes

**Missing\_geno**: check for missingness in SNP calling, remove SNPs with a low genotyping rate.

**genotype\_missingness\_check**: _Not sure what this does, i think it links with the above but could be an individual missingness check_

**MAF\_summary**: Do you want to check the minor allele frequency in the summary stats dataset?

**MAF\_threshold**: threshold of SNPs to remove based on their minor allele frequency. Every SNP below this threshold will be removed

**MAF\_genotype**: Do you want to check the minor allele frequncy in the genotype test dataset?

**INFO\_summary**: INFO score QC?

**INFO\_threshold**: What INFO score threshold do you want to use? Every SNP below this threshold will be remove

**SE\_summary**: do you want to QC on standard error? 

**SE\_threshold**: what SE threshold do you want to use? Every SNP above this threhold will be removed.

**Raven\_out\_info\_directory**: DEFUNCT?

**Chromosomes\_to\_analyse**: which chromosomes do you want to analyse? Must be in the form of an array.

**p1**: p-value threshold one for intelligent pruning in plink

**p2**: p-value threshold two for intelligent pruning in plink

**r2**: r2 value threhold for intelligent pruning in plink

**window**: size of window for intelligent pruning in plink. scale to KB.

**sig\_thresholds**: p-value thresholds for Polygenic risk score. must be an array.

**sig\_thresholds\_lower\_bounds**: element-by-element lower bound for match-up to plink's p value threshold value. must equal the length of sig\_thresholds and also be an array.

**Extra\_analysis**: Do you want to perform more than the standard whole genome polygenic risk score?

**Name\_of\_extra\_analysis**: State either "Genes", "Pathways" or an array of (Genes Pathways)

**randomise**: do you want to perform permutation analysis on a pathway-specific polygenic risk score?

**permutations**: number of iterations to permute the pathway specific polygenic risk score.

**Magma\_validation\_set\_name**: UNSURE

**Gene\_regions**: use gene boundaries as limits or include extended transcription regions? either normal, extended or both.

**whole\_genome\_genic**: perform the whole genome polygenic risk score wihthin genic limits?

**Gene\_specific\_PRS**: perform a polygenic risk score per gene?

