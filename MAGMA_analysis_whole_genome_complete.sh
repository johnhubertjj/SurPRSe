#! /bin/bash

whereami=$(uname -n)
echo "$whereami"

if [[ "$whereami" == *"raven"* ]]; then

#Load Magma
module load magma/1.06
fi

if [[ "${whereami}" == 'v1711-0ab8c3db.mobile.cf.ac.uk' ]]; then
cd ~/Dropbox/whole_genome_testing/output/
echo pwd
echo "i'm here"
fi

# Read in the significance thresholds you wish to perform on the dataset: make sure you read the variable in as an array 
# Gets a bit difficult with arrays in bash, the argument could be any length and you need to print each element in as an argument
# so assign sig to all the arguments to the script
sig=( "$@" )
echo "${sig[@]}"

# take the length of the array
length_of_total_array=${#sig[@]}

# you know the length of the other arguments but is also zero based so add an argument = length_of_arguments + 1 so that you can loop through the significance thresholds effectively
num6=6
num5=5
num4=4
num3=3
num2=2
num1=1

length_of_sig=$((${length_of_total_array} - ${num5}))

# perform a sum of the arguments to exclude the significance threshold arguments
length_of_total_array_for_validation=${length_of_sig}
length_of_total_array_for_training=$((${length_of_sig} + ${num1}))
length_of_total_array_for_validation_set_name_MAGMA=$((${length_of_sig} + ${num2}))
length_of_total_array_for_whether_to_perform_magma=$((${length_of_sig} + ${num3}))
length_of_total_array_for_Gene_regions=$((${length_of_sig} + ${num4}))

# assign the other arguments here
validation_set_name=${sig[@]:${length_of_total_array_for_validation}:${num1}}
training_set_name=${sig[@]:${length_of_total_array_for_training}:${num1}}
validation_set_name_MAGMA=${sig[@]:${length_of_total_array_for_validation_set_name_MAGMA}:${num1}}
Perform_magma_analysis=${sig[@]:${length_of_total_array_for_whether_to_perform_magma}:${num1}}
Gene_regions=${sig[@]:${length_of_total_array_for_Gene_regions}:${num1}}

# Delete the other arguments from the significance threshold array argument
delete=( ${validation_set_name} ${training_set_name} ${validation_set_name_MAGMA} ${Perform_magma_analysis} ${Gene_regions} )

for del in ${delete[@]}
do
	sig=("${sig[@]/$del}")
done
echo ${sig[@]}

# Create a variable that will allow you to loop through the array correctly as it is zero indexed 
length_of_thresholds=`echo "${#sig[@]}"` 
echo ${length_of_thresholds}
length_of_thresholds=`echo "$((${length_of_thresholds} - $num6))"`
echo ${length_of_thresholds}

# Whether you wish to include regulatory regions or not

if [[ ${Gene_regions} == normal || ${Gene_regions} == both ]]; then
echo "works"
	if [ ${Perform_magma_analysis} == TRUE ]; then
echo "works2"
		for i in `seq 0 ${length_of_thresholds}` ;
		do 
			echo ./MAGMA_set_analysis/${sig[i]}${validation_set_name}_${training_set_name}_P_Vals_reference_file_magma_analysis.txt
			echo ${sig[@]}
			echo ${length_of_thresholds}
			awk '{ print $2 }' ./MAGMA_set_analysis/${sig[i]}${validation_set_name}_${training_set_name}_P_Vals_reference_file_magma_analysis.txt > ./MAGMA_set_analysis/${sig[i]}${validation_set_name}_${training_set_name}_Filter_SNPS_per_pvalue_magma_analysis.txt
			magma --annotate filter=./MAGMA_set_analysis/${sig[i]}${validation_set_name}_${training_set_name}_Filter_SNPS_per_pvalue_magma_analysis.txt --snp-loc ${validation_set_name_MAGMA} --gene-loc NCBI37.3.gene.loc --out ./MAGMA_set_analysis/${sig[i]}_${validation_set_name}_${training_set_name}_SNPs_whole_genome_magma
			magma --bfile ${validation_set_name_MAGMA} --gene-annot ./MAGMA_set_analysis/${sig[i]}_${validation_set_name}_${training_set_name}_SNPs_whole_genome_magma.genes.annot --covar file=CLOZUK2.r7.select2PC.eigenvec.txt include=PC1-PC5,PC6,PC9,PC11,PC12,PC13,PC19 --out ./MAGMA_set_analysis/${sig[i]}gene_annotation_for_CLOZUK_without_selecting_genes_magma
		done 
	else
		magma --annotate --snp-loc ${validation_set_name}_${training_set_name}_FULL_GENOME.bim --gene-loc NCBI37.3.gene.loc --out ${validation_set_name}_${training_set_name}_SNPs_whole_genome
	fi	
fi
 
if [[ ${Gene_regions} == "extended" || ${Gene_regions} == "both" ]]; then

	if [[ ${Perform_magma_analysis} == "TRUE" ]]; then

		for i in `seq 0 ${length_of_thresholds}` ;
		do 
			awk '{print $2}' ./MAGMA_set_analysis/${sig[i]}${validation_set_name}_${training_set_name}_P_Vals_reference_file_magma_analysis.txt > ./MAGMA_set_analysis/${sig[i]}${validation_set_name}_${training_set_name}_Filter_SNPS_per_pvalue_magma_analysis.txt
			magma --annotate window=35,10 filter=./MAGMA_set_analysis/${sig[i]}${validation_set_name}_${training_set_name}_Filter_SNPS_per_pvalue_magma_analysis.txt --snp-loc ${validation_set_name_MAGMA} --gene-loc NCBI37.3.gene.loc --out ./MAGMA_set_analysis/${sig[i]}_${validation_set_name}_${training_set_name}_SNPs_whole_genome_magma_extended
			magma --bfile ${validation_set_name_MAGMA} --gene-annot ./MAGMA_set_analysis/${sig[i]}_${validation_set_name}_${training_set_name}_SNPs_whole_genome_magma_extended.genes.annot --covar file=CLOZUK2.r7.select2PC.eigenvec.txt include=PC1-PC5,PC6,PC9,PC11,PC12,PC13,PC19 --out ./MAGMA_set_analysis/${sig[i]}gene_annotation_for_CLOZUK_without_selecting_genes_magma
		done 
	else
		magma --annotate window=35,10 --snp-loc ${validation_set_name}_${training_set_name}_FULL_GENOME.bim --gene-loc NCBI37.3.gene.loc --out ${validation_set_name}_${training_set_name}_SNPs_whole_genome_extended
	fi
fi	
