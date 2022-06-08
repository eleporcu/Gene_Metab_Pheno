# Gene_Metab_Pheno
Codes used in "Exploiting the mediating role of the metabolome to unravel transcript-to-phenotype associations" (Auwerx et al, bioRxiv 2022)

######## MR_uni.R

MR_uni.R requires as input a matrix containing the following columns: 
-SNPs: SNPs to be included in the model; 
-GENE: the univariate effect size of these SNPs on gene expression (these estimates come from an eQTL study) ; 
-GENE_N: the sample size of the eQTL study; 
-BETA_GWAS: the univariate effect sizes of these SNPs on the phenotype; 
-N: the sample size of the GWAS.

How to run the script:
R < MR_uni.R --no-save ENSG00000115216_matrix.txt

The output is a file containing the following columns: 
-gene: name of the gene tested;
-alpha_ORIGINAL: causal effect estimated before applying heterogeneity test;
-SE_ORIGINAL: standard error of alpha_ORIGINAL; 
-P_ORIGINAL: P-value calculated from alpha_ORIGINAL and SE_ORIGINAL; 
-N_ORIGINAL: number of SNPs used as instruments; 
-P_het_ORIGINAL: P-value of heterogenity test; 
-alpha: causal effect estimated after removing SNPs detected as outliers by the heterogenity test; 
-SE: standard error of alpha; 
-P: P-value calculated from alpha and SE; 
-N: number of SNPs left after outlier removal; 
-Phet: P-value of heterogenity test after outlier removal; 
-N_outlier: number of removed outliers.


######## MR_multi.R

MR_multi.R requires as input a matrix containing the following columns: 
-SNPs: SNPs to be included in the model; 
-GENE: the univariate effect size of these SNPs on gene expression (these estimates come from an eQTL study); 
-METABOLITE: the univariate effect size of these SNPs on the metabolite (these estimates come from an mQTL study);
 -BETA_GWAS: the univariate effect sizes of these SNPs on the phenotype.

How to run the script:
R < MR_multi.R --no-save ENSG00000261015_M32654_bilirubin_matrix.txt

The output gives the direct effect


######## power.R 

Code we used for the simulation analysis performed to estimate differences in power when ussing the mediation analysis, as opposed to transcriptome-wide Mendelian rrandomization (TWMR).
How to run the script: R < power.R --no-save


How to run the script:
R < power.R --no-save

