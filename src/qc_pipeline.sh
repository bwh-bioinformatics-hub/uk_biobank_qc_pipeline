#!bin/bash

## step1: extract subjects with IDP (n~=40k) and gender matched
# prepare the idp_iid.fam
# TOADD: code to generate idp_iid.fam

# this step is common for all chr, so we just need to do once here. 
[ -f idp_iid.genderMatched.fam ] || {
        # get male from genotype 
        plink2 --bgen genotype/ukb_imp_chr${1}_v3.bgen ref-first \
         --keep-males \
         --sample genotype/ukb28142_imp_chr22_v3_s487296.sample \
         --make-just-fam \
         --out qc_output/imp_chr22_step4.male
        # male in both genotype and phenotype
        sort <(cut -f1 qc_output/imp_chr22_step4.male.fam) <(awk '$4==1{print $2}' genotype/ukb28142_imp_chr22_v3_s*.sample) | uniq -d > qc_output/validated_male.id

        # do the same for female
	plink2 --bgen genotype/ukb_imp_chr${1}_v3.bgen ref-first \
         --keep-females \
         --sample genotype/ukb28142_imp_chr22_v3_s487296.sample \
         --make-just-fam \
         --out qc_output/imp_chr22_step4.female
        # male in both genotype and phenotype
        sort <(cut -f1 qc_output/imp_chr22_step4.female.fam) <(awk '$4==2{print $2}' genotype/ukb28142_imp_chr22_v3_s*.sample) | uniq -d > qc_output/validated_female.id
        
	# extract only the gender-matched subjects
        cat qc_output/validated_male.id qc_output/validated_female.id | fgrep -w -f - idp_iid.fam > idp_iid.genderMatched.fam
	
	awk '{print $1 " " $1}' idp_iid.genderMatched.fam > idp_iid.genderMatched1.fam
	rm idp_iid.genderMatched.fam
	mv idp_iid.genderMatched1.fam idp_iid.genderMatched.fam
}

## step2: remove subjects with >1% genotype missingness

plink2 \
 --bgen genotype/ukb_imp_chr${1}_v3.bgen ref-first \
 --sample genotype/ukb*_imp_chr${1}_*.sample \
 --keep idp_iid.genderMatched.fam \
 --mind 0.01 \
 --export bgen-1.2 \
 --make-just-fam \
 --out qc_output/imp_chr${1}_step2

#generate sample file for gender matched samples
fgrep -wf <(cut -f1 qc_output/imp_chr${1}_step2.fam) genotype/ukb28142_imp_chr${1}_v3_*.sample > qc_output/imp_chr${1}.sample
sed -i '1s/^/ID_1 ID_2 missing sex\n0 0 0 D\n/' qc_output/imp_chr${1}.sample

## step3: only keep bi-allelic SNPs
## step4: remove SNPs with low genotyping rate (--geno 0.01)
## step5: remove SNPs with HWE ... ()
## step6: MAF
## step7: Perform LD-based pruning to remove highly-correlated SNPs

plink2 \
 --bgen qc_output/imp_chr${1}_step2.bgen ref-last \
 --sample qc_output/imp_chr${1}.sample \
 --snps-only --min-alleles 2 --max-alleles 2 \
 --rm-dup exclude-all \
 --geno 0.01 \
 --hwe 1e-6 \
 --maf 0.01 \
 --write-snplist \
 --indep-pairwise 200 50 0.25 \
 --out qc_output/imp_chr${1}_step7

# debug
echo "SNP number in write-snplist:" `wc -l qc_output/imp_chr${1}_step7.snplist`
echo "SNP number in pruned lists:" `wc -l qc_output/imp_chr${1}_step7.prune.*`

## step8: Heterozygosity check

plink2 \
 --bgen qc_output/imp_chr${1}_step2.bgen ref-last \
 --sample qc_output/imp_chr${1}.sample \
 --extract qc_output/imp_chr${1}_step7.prune.in \
 --het \
 --out qc_output/imp_chr${1}_step8

# valid samples: F-coefficient is within +/-3 sd of the overall mean F-coefficient.
# output: qc_output/imp_chr${1}_step8.valid.sample  
Rscript src/het.r qc_output/imp_chr${1}_step8.het qc_output/imp_chr${1}_step8.valid.sample

## step9: check sample kinship using KING-robust (similar to IBS/IBD check)
## step10: run population stratification using pca
plink2 \
 --bgen qc_output/imp_chr${1}_step2.bgen ref-last \
 --sample qc_output/imp_chr${1}.sample \
 --extract qc_output/imp_chr${1}_step7.prune.in \
 --keep qc_output/imp_chr${1}_step8.valid.sample \
 --king-cutoff 0.125 \
 --pca \
 --out qc_output/imp_chr${1}_step10

## step10: check mismatching SNPs btw base (GWAS) and target (bgen) data
# output: 
Rscript src/mismatch.r ${1}

## step11: extract the final set of SNPs and samples, write to the final bgen file

plink2 \
 --bgen qc_output/imp_chr${1}_step2.bgen ref-last \
 --sample qc_output/imp_chr${1}.sample \
 --exclude qc_output/imp_chr${1}.mismatch \
 --extract qc_output/imp_chr${1}_step7.snplist \
 --keep qc_output/imp_chr${1}_step10.king.cutoff.in.id \
 --out qc_output/imp_chr${1}.final \
 --export bgen-1.2


