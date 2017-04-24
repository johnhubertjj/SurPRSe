#! /Library/Frameworks/Python.framework/Versions/2.7/bin/python2.7

# -------------- READ IN FUNCTIONS -------------- #

import pandas as pd
import numpy as np
from decimal import Decimal
from multiprocessing import Pool
import time
import os
import platform
import argparse
import fileinput
import csv

# set the script for analysis depending on the system on which you are running the PRS analysis
platform_system = platform.uname()

if platform_system[1] == 'v1711-0ab8c3db.mobile.cf.ac.uk' :
    os.chdir('/Users/johnhubert/Dropbox/whole_genome_testing/output/')
    chromosome_number = 22
elif platform_system[1] == 'JJ' :
    os.chdir('~/Documents/testing_/output')
    chromosome_number = 22
elif platform_system[1] == 'rocks.psycm.cf.ac.uk' :
    os.chdir('~/output/')
    # below is wrong, raven will not always be on raven13 #
    chromosome_number = os.getenv('PBS_ARRAY_INDEX')
elif platform_system[1] == 'raven13' :
    os.chdir('/scratch/c1020109/PR54/PGC_CLOZUK_PRS/output/')
    chromosome_number = os.getenv('PBS_ARRAY_INDEX')
#
os.getcwd()

timer_start = time.time()

# -------------- READ IN DATA -------------- #

# Read in Annot Data #
annot_names = ["CHR", "SNP", "BP", "GENE", "BP_START", "BP_END"]
annot = pd.read_table('MAGMA_Gene_regions_for_python_script.txt', sep=' ', names=annot_names)

# Read in CLOZUK_PGC Assoc Data Pruned Data # #Is this the PGC alternative?#

CLOZUK_PGC_names = ["ROW.NUM", "ROW.ID", "CHR", "SNP", "BP", "A1", "A2", "BETA", "P", "MAF"]
CLOZUK_PGC_clump = pd.read_table('CLOZUK_GWAS_BGE_CLUMPED_PGC_MAF_FINAL.txt', sep=' ', names=CLOZUK_PGC_names)

# If the data is from R, the indexing standard needs to be changed
CLOZUK_PGC_clump[['ROW.NUM']] = CLOZUK_PGC_clump[["ROW.NUM"]] - 1

# File that Emily made herself: took the MAF from GERAD and then applied it above.
# For me I would take the Minor allele frequency from CLOZUK and then apply the BETA and p for PGC for the common SNPs