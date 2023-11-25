plink --vcf genotypes.vcf --pca 10

plink --vcf genotypes.vcf --freq

plink --vcf genotypes.vcf --linear --pheno CB1908_IC50.txt --covar plink.eigenvec --allow-no-sex --out phenotype_gwas_results