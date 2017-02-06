#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python2.7

# -------------- READ IN FUNCTIONS -------------- #

import pandas as pd
import numpy as np
from decimal import Decimal
from multiprocessing import Pool
import time
import os
import platform
import fileinput

# -------------- DECLARE VARIABLES -------------- #

indiv_snp = []
significance_threshold = os.getenv('Significance_threshold')

# set the script for analysis depending on the system on which you are running the PRS analysis
platform_system = platform.uname()

if platform_system[1] == 'v1711-0ab8c3db.mobile.cf.ac.uk' :
    os.chdir('/Users/johnhubert/Documents/testing_PRS_chromosome_22/output/')
    chromosome_number = 22
elif platform_system[1] == 'JJ' :
    os.chdir('~/Documents/testing_PRS_chromosome_22/output')
    chromosome_number = 22
elif platform_system[1] == 'rocks.psycm.cf.ac.uk' :
    os.chdir('~/output/')
    # below is wrong, haven't run a loop on rocks yet #
    chromosome_number = os.getenv('PBS_ARRAY_INDEX')
elif platform_system[1] == 'raven13' :
    os.chdir('/scratch/c1020109/PR54/PGC_CLOZUK_PRS/output/')
    chromosome_number = os.getenv('PBS_ARRAY_INDEX')

os.getcwd()

timer_start = time.time()

# --------------- CREATE SIGNIFICANCE THRESHOLD FUNCTION --------------- #
def Sig_func() :
    sig_threshold = np.array(sys.argv)
    to_exclude = [0]
    sig_threshold_without_scriptpath = np.delete(sig_threshold, to_exclude)
    print(sig_threshold_without_scriptpath)

    for i in range(0, len(sig_threshold_without_scriptpath)) :
    # alter file input here for each and assign to a new object
    # make sure to keep is separate from other sources of code as well


# -------------- CREATE PRS FUNCTION -------------- #

# Define function #
def PRS_func(gene):

    # Create empty datasets for results #
    indiv_id = range(1, 35802)
    results=pd.DataFrame(index = indiv_id, columns = ['GeneName'])

    # Find SNPs in gene #
    gene_snp = annot[annot.GENE == unigene.iloc[gene]['GENE']]
    print(gene_snp)
    # Delete repeat data #
    del gene_snp['CHR']
    del gene_snp['BP']

    # Find SNP nos and effect sizes corresponding to SNPs in genes- also finds pruned SNPs per gene #
    effect = pd.merge(gene_snp, CLOZUK_PGC_clump, on = ['SNP'])
    print(effect)

    # Read in Genetic Data per Subject Pruned Data #
    if (len(effect.index) != 0):
        # The raw file uses recode-A NOT AD and the ROW ID is a recorder of the row for the right SNP (+6 for python because of zero indexing but +5 otherwise)
        # then take the headers off - command line will work on shell script
        gen_clump = open('CLOZUK_GWAS_BGE_CLUMPED_chr' + str(chromosome_number) + '_no_head.raw', 'r')

        indiv = 0
        
        for line in gen_clump:
            
            # Individual number for results table #
            indiv = indiv + 1
        
            # Create vector for individual genotype values #
            indiv_snp = line.split()
            
            # Find PRS score for each individual #

            prs=0.00
            
            # Loop through effect dataset #
            for i in range(0, len(effect.index)):
                # Careful here, in R the indexing starts at 1 so the rows will be out out order by one
                if indiv_snp[(effect.iloc[i]['ROW.NUM']+6)] != 'NA':
                    # Calculate the PRS, find individual genotype corresponding to correct SNP using rowID #
                    prs = prs + (float(indiv_snp[(effect.iloc[i]['ROW.NUM'] + 6)]) * effect.iloc[i]['BETA'])
                else:
                    prs = prs + (2*effect.iloc[i]['BETA'] * effect.iloc[i]['MAF'])
        
                results.loc[indiv, 'GeneName'] = (prs / float(len(effect.index)))

        gen_clump.close()

        # Rename the column header with the actual gene name #
        results = results.rename(columns = {'GeneName': unigene.iloc[gene]['GENE']})

        # Send results to main processor #
        return results

# If on the main processor #
if __name__ == '__main__':

    # -------------- READ IN DATA -------------- #

    # Read in Unique Gene file #
    unigene_names = ["GENE"]
    unigene = pd.read_table('PGC_CLOZUK_unique_genes.txt', names = unigene_names)

    # Read in Pruned Bim file #
    bim_names = ["CHR", "SNP", "GD", "BP", "A1", "A2"]
    bim_clump = pd.read_table( 'CLOZUK_GWAS_BGE_CLUMPED_chr' + str(chromosome_number) + '.bim', names=bim_names)

    
    # Read in Annot Data #
    annot_names = ["CHR", "SNP", "BP", "GENE", "BP_START","BP_END"]
    annot = pd.read_table('/output/MAGMA_Gene_regions_for_python_script_chr_' + str(chromosome_number) + '.txt', sep=' ', names=annot_names)
    
    # Read in CLOZUK_PGC Assoc Data Pruned Data # #Is this the PGC alternative?#
    
    CLOZUK_PGC_names = ["ROW.NUM","ROW.ID", "CHR", "SNP", "BP", "A1", "A2", "BETA", "P", "MAF"]
    CLOZUK_PGC_clump = pd.read_table('./output/CLOZUK_GWAS_BGE_CLUMPED_PGC_MAF_FINAL' + str(chromosome_number) + '.txt', sep = ' ', names = CLOZUK_PGC_names)

    # If the data is from R, the indexing standard needs to be changed
    CLOZUK_PGC_clump[['ROW.NUM']] = CLOZUK_PGC_clump[["ROW.NUM"]] - 1

    # File that Emily made herself: took the MAF from GERAD and then applied it above.
    # For me I would take the Minor allele frequency from CLOZUK and then apply the BETA and p for PGC for the common SNPs

    # -------------- CALL FUNCTION IN PARALLEL -------------- #

    pool = Pool()
    res = pool.map(PRS_func, range(10))

    # -------------- GROUP AND OUTPUT RESULTS -------------- #

    prs_results_clump = pd.concat(res, axis = 1)

    prs_results_clump.to_csv('./PRS_scoring/prs_results_clump_miss' + str(chromosome_number), header = True, index = None, sep = ' ')

timer_end=time.time()-timer_start

print timer_end

# At the moment the script loops through quite a few things and requires the input of a large dataset (albeit line by line) into the python environment
    # Improvements: make it work for the plink format: read in and parse through the binary file so that it becomes more applicable and easier to run
    # will save on memory, also try not to loop so many times, it's effectively adding the PRS each time and then overwriting that vector, would be better to pre-create the PRS vector first and then add up at the end in one function instead of overwriting
    # It does however, do this in parallel, (which is good) but it will need to be tracked so that it can work