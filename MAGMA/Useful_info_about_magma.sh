#! /bin/bash

sig=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)
length_PGC=(521 1435 4787 13603 23369 40794 57657 74106 89961)

#!/bin/bash

cd ~/Dropbox/testing_PRS_chromosome_22/
chromosome_number=22

sig=(0.0001 0.001 0.01 0.05 0.1 0.2 0.3 0.4 0.5)

for i in `seq 0 8` ;
do

#wc -l ${sig[i]}update_Bed_files_for_MAGMA_analysis.txt >> Number_of_SNPs_before_MAGMA.txt
wc -l ${sig[i]}CLOZUK_PGC_P_Vals_reference_file_using_PRS.txt >> Number_of_SNPs_before_MAGMA_2.txt
#wc -l ${sig[i]}unique_merged_PRS_MAGMA_COMMON_GENES.txt >> Number_of_Genes_before_MAGMA.txt
done

