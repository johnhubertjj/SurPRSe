#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python2.7

# -------------- DECLARE VARIABLES -------------- #

indiv_snp=[]

# -------------- READ IN FUNCTIONS -------------- #

import pandas as pd
import numpy as np
from decimal import Decimal
from multiprocessing import Pool
import time

timer_start=time.time()


# -------------- CREATE PRS FUNCTION -------------- #

# Define function #
def f(gene):

    # Create empty datasets for results #
    indiv_id=range(1, 13165)
    results=pd.DataFrame(index=indiv_id, columns=['GeneName'])

    # Find SNPs in gene #
    gene_snp=annot[annot.GENE==unigene.iloc[gene]['GENE']]
    print(gene_snp)
    # Delete repeat data #
    del gene_snp['CHR']
    del gene_snp['BP']

    # Find SNP nos and effect sizes corresponding to SNPs in genes- also finds pruned SNPs per gene #
    effect=pd.merge(gene_snp, igap_clump, on=['SNP'])
    print(effect)

    # Read in Genetic Data per Subject Pruned Data #
    if (len(effect.index)!=0):
        # The raw file uses recode-A NOT AD and the ROW ID is a recorder of the row for the right SNP (+6 for python because of zero indexing but +5 otherwise)
        # then take the headers off - command line will work on shell script
        gen_clump = open('Clumping/GERAD_geno_GeneSNPs_igapnogerad_clump_nohead.raw', 'r')

        indiv=0
        
        for line in gen_clump:
            
            # Individual number for results table #
            indiv = indiv+1
        
            # Create vector for individual genotype values #
            indiv_snp = line.split()
            
            # Find PRS score for each individual #

            prs=0.00
            
            # Loop through effect dataset #
            for i in range(0, len(effect.index)):
                
                if indiv_snp[(effect.iloc[i]['ROW_ID']+6)]!= 'NA':
                    # Calculate the PRS, find individual genotype corresponding to correct SNP using rowID #
                    prs = prs+(float(indiv_snp[(effect.iloc[i]['ROW_ID']+6)])*effect.iloc[i]['BETA'])
                else:
                    prs = prs+(2*effect.iloc[i]['BETA']*effect.iloc[i]['MAF'])
        
                results.loc[indiv, 'GeneName'] = (prs / float(len(effect.index)))

        gen_clump.close()

        # Rename the column header with the actual gene name #
        results=results.rename(columns = {'GeneName': unigene.iloc[gene]['GENE']})

        # Send results to main processor #
        return results

# If on the main processor #
if __name__=='__main__':

    # -------------- READ IN DATA -------------- #

    # Read in Unique Gene file #
    unigene_names=["GENE"]
    unigene=pd.read_table('PGC_CLOZUK_unique_genes.txt', names=unigene_names)

    # Read in Pruned Bim file #
    bim_names=["CHR", "SNP", "GD", "BP", "A1", "A2"]
    bim_clump=pd.read_table('CLOZUK_GWAS_BGE_CLUMPED_chr22.bim', names=bim_names)

    
    # Read in Annot Data #
    annot_names=["CHR", "SNP", "BP", "GENE", "BP_START","BP_END"]
    annot=pd.read_table('all_genes_with_SNPS.txt', sep=' ', names=annot_names)
    
    # Read in IGAP_noGERAD Assoc Data Pruned Data # #Is this the PGC alternative?#
    
    igap_names=["ROW.NUM","ROW.ID", "CHR", "SNP", "BP", "A1", "A2", "BETA", "P", "MAF"]
    igap_clump=pd.read_table('CLOZUK_GWAS_BGE_CLUMPED_PGC_MAF_FINAL.txt', sep=' ', names=igap_names)
    # File that Emily made herself: took the MAF from GERAD and then applied it above.
    # For me I would take the Minor allele frequency from CLOZUK and then apply the BETA and p for PGC for the common SNPs

    # -------------- CALL FUNCTION IN PARALLEL -------------- #

    pool=Pool()
    res=pool.map(f, range(10))

    # -------------- GROUP AND OUTPUT RESULTS -------------- #

    prs_results_clump=pd.concat(res, axis=1)

    prs_results_clump.to_csv('PRS_scoring/prs_results_clump_miss', header=True, index=None, sep=' ')

timer_end= time.time()-timer_start

print timer_end


