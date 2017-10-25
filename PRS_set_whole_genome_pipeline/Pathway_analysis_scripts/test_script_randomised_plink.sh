#! /bin/bash

pathways=(5HT_2C abnormal_behavior Cav2_channels FMRP_targets Lek2015_LoFintolerant_90 abnormal_grooming_behavior abnormal_long_term_potentiation abnormal_nervous_system_electrophysiology Calcium_ion_import_GO0070509 Membrane_depolarization_during_action_potential_GO0086010 Synaptic_transmission_GO0007268) 

sig_thresholds=(0.05 0.5 1)


Pathway_output_directory="/home/johnhubert/Documents/testing_random_gene_sets/Profiles/"
Gene_output_directory="/home/johnhubert/Documents/testing_random_gene_sets/Scores/"
bim_file_output_directory="/home/johnhubert/Documents/ALSPAC_gene_pathway_pipeline_test/CLOZUK_PGC2noclo_ALSPAC_output/Genes/"

# Create PRS profiles for each significance threshold specified

for setname in "${pathways[@]}" ;
do

	echo ${setname}
	
for i in "${sig_thresholds[@]}" ;
do

for l in `seq 1 1000` ;
do
	echo ${l}

	filename="${Pathway_output_directory}${setname}_random_${l}_with_${i}"
	echo ${filename}

	plink --bfile ${bim_file_output_directory}ALSPAC_CLOZUK_PGC2noclo_normal_gene_regions_Clumped_whole_genome_final --score ${Gene_output_directory}${setname}_random_${l}_with_${i}.score --out ${filename}


done
done
done


