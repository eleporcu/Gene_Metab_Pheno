# Gene_Metab_Pheno
Codes used in "Exploiting the mediating role of the metabolome to unravel transcript-to-phenotype associations" (Auwerx et al, bioRxiv 2022)


MR_uni.R requires as input a matrix containing the univariate effect size of n SNPs on gene expression (these estimates come from an eQTL study) 
and the univariate effect sizes on the phenotype. The last three columns are: BETA_GWAS SE N

How to run the script:

R < MR_uni.R --no-save bmi.matrix.betaGWAS bmi genes.N

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



