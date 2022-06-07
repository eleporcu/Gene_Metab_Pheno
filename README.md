# Gene_Metab_Pheno
Codes used in "Exploiting the mediating role of the metabolome to unravel transcript-to-phenotype associations" (Auwerx et al, bioRxiv 2022)

########   MR_uni.R
MR_uni.R requires as input a matrix containing as columns:
-SNPs: the SNPs to be included in the model
-GENE: the univariate effect size of the SNPs on gene expression (these estimates come from an eQTL study) 
-GENE_N: the sample size of the eQTL study
-BETA_GWAS: the univariate effect sizes of the SNPs on the phenotype 
-N: the sample size of the GWAS

How to run the script:

R < MR_uni.R --no-save ENSG00000115216_matrix.txt

The output files is a file containing the following columns: 
-gene: name of the gene tested 
-alpha_ORIGINAL: causal effect estimated by revTWMR before applying heterogeneity test 
-SE_ORIGINAL: standard error of alpha_original 
-P_ORIGINAL: Pvalue calculated from alpha_original and SE_original 
-N_ORIGINAL: number of SNPs used as instruments 
-P_het_ORIGINAL: Pvalue of heterogenity test 
-alpha: causal effect estimated after removing the SNPs detected as outliers by the heterogenity test 
-SE: standard error of alpha 
-P: Pvalue calculated from alpha and SE 
-N: number of SNPs left after outlier removal 
-Phet: Pvalue of heterogenity test after outlier removal
-N_outlier: number of removed outliers





