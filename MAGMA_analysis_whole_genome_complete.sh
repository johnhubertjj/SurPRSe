#! /bin/bash
module load magma/versionhere

# Read in the significance thresholds you wish to perform on the dataset: make sure you read the variable in as an array 
sig=( "${1[@]}" )

# Convert arguments into human readable variables
validation_set_name=$2
training_set_name=$3
Validation_dataset=$4
Perform_magma_analysis=$5

# Create a variable that will allow you to loop through the array correctly as it is zero indexed 
length_of_thresholds=`echo "${#sig[@]}"` 
num1=1
length_of_thresholds=`echo "$((${length_of_array} - ${num1}))"`


if [ ${Perform_magma_analysis} == "TRUE"]; then

	for i in `seq 0 ${length_of_array}` ;
	do 
		awk '{print $2}' ./MAGMA_set_analysis/${sig[i]}${validation_set_name}_${training_set_name}_P_Vals_reference_file_magma_analysis.txt > ./MAGMA_set_analysis/${sig[i]}${validation_set_name}_${training_set_name}_Filter_SNPS_per_pvalue_magma_analysis.txt
		magma --annotate filter=./MAGMA_set_analysis/${sig[i]}${validation_set_name}_${training_set_name}_Filter_SNPS_per_pvalue_magma_analysis.txt --snp-loc ${validation_set_name}_${training_set_name}_MAGMA_FULL_GENOME_.bim --gene-loc NCBI37.3.gene.loc --out ./MAGMA_set_analysis/${sig[i]}_${validation_set_name}_${training_set_name}_SNPs_whole_genome_magma
		magma --bfile ${validation_set_name}_${training_set_name}_MAGMA_FULL_GENOME --gene-annot ./MAGMA_set_analysis/${sig[i]}_${validation_set_name}_${training_set_name}_SNPs_whole_genome_magma.genes.annot --covar file=CLOZUK2.r7.select2PC.eigenvec.txt include=PC1-PC5,PC6,PC9,PC11,PC12,PC13,PC19 --out ./MAGMA_set_analysis/${sig[i]}gene_annotation_for_CLOZUK_without_selecting_genes_magma
	done 
else
	for i in `seq 0 ${length_of_array}` ;
	do 
		awk '{print $2}' ${sig[i]}${validation_set_name}_${training_set_name}_P_Vals_reference_file.txt > ${sig[i]}${validation_set}_${training_set_name}_Filter_SNPS_per_pvalue.txt
		magma --annotate filter=${sig[i]}CLOZUK_PGC_Filter_SNPS_per_pvalue.txt --snp-loc CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS.bim --gene-loc NCBI37.3.gene.loc --out ${sig[i]}_CLOZUK_PGC_SNPs_whole_genome
		magma --bfile CLOZUK_PGC_FULL_GENOME_without_erroneous_SNPS --gene-annot ${sig[i]}_CLOZUK_PGC_SNPs_whole_genome.genes.annot --covar file=CLOZUK2.r7.select2PC.eigenvec.txt include=PC1-PC5,PC6,PC9,PC11,PC12,PC13,PC19 --out ${sig[i]}gene_annotation_for_CLOZUK_without_selecting_genes
	done 
fi	

