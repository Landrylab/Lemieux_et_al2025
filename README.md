# Lemieux et al. 2025 code and Supplementary Material

The repository is divided in the following 5 directories

## /cytometry
Contains the R script used to analyse the flow cytometry experiments and the compressed files for both cytometry experiments.
	analysis_cytometry.R \
	GFP_XX.zip

## /growthcurves_analysis
Contains the R scripts to analyse the 4 PCA experiments 
	AUC_analysis_XX.R

### /gc_data
Contains the area under the curve and derivative growth rate computed, one for each PCA experiments. \
	condensed_XX.csv \
There is also a file listing the affinity data retrieved from previous study used for comparison with PCA data. \
	ref_affinity.xlsx

## /PBD_selection
Contains the R scripts and jupyter notebook used for the selection of the PBD and peptides from the PRM DB. The files are numbered in the order of execution to obtain the final sequences used in the last DHFR PCA experiment. 

### /reference_data
Contains the files used as input in the script of the parent directory. There are 2 files used for the codon optimization for yeast of dna sequences, dict_condon.csv and dict_usage.csv, used in 4_Codon_optimization.ipynb. \
There are also 2 files, PRM_Master.csv and PWM_updated.rds, extracted from the PRM-DB. Here is their documentation : https://prm-db.org/download.php

## /supplementary_material
Contains the Supplementary File 1 - detailed methods - and Supplementary File 2 - supplementary Tables 1 to 6. \
	Supplementary_File1.pdf \
	Supplementary_File2.xlsx

## WB_images.zip
Compressed folder with the original images from the Western Blots and Ponceau colorations.
	
See Lemieux et al. 2025 for further details. 
