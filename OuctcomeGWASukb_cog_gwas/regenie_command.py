
import os
import pandas as pd

os.system("module load PLINK/2.00a2.3_x86_64")
os.system("module load Regenie/3.2.4-GCC-11.2.0")
os.system("module load bgen/1.1.6")
os.system("module load BGEN-enkre/1.1.7-GCC-11.2.0")


cog_phenotype="/edgehpc/dept/compbio/human_genetics/users/cchen11/ukbb_cogexome/data/pheno_v3/ukbb500k.fiall_baseline_model2.EUR.regenie.v3.txt"
UKB_PPP_fam="/edgehpc/dept/compbio/human_genetics/users/bsun/UKB_PPP_pQTLs/analysis/interim/ukb_genetic_variants/ukb_imp_chr1_mac50_info07_b0_7_patched_bfiles.fam" 
geenotype_file="/edgehpc/dept/compbio/human_genetics/ukbb_temp/Release2/Genotyped/ukb_chrauto_v2_snpqc"



##To find the family and individual ids of the file ukbb500k
os.system("plink2 --pfile  /edgehpc/dept/compbio/human_genetics/ukbb_temp/Release2/Genotyped/ukb_chrauto_v2_snpqc --make-bed --out ukb_chrauto_v2_snpqc")


#identify the common samples present in ukbb500k and UKB_PPP
UKB_PPP_fam_df=pd.read_csv(UKB_PPP_fam,sep="\t",header=None)
UKB_PPP_fam_df.columns=["ukpp_FID", "ukpp_IID","ukpp_PID","ukpp_MID","ukpp_sex","ukpp_disease"]

ukbb500k_fam_df=pd.read_csv("ukb_chrauto_v2_snpqc.fam",sep="\t",header=None)
ukbb500k_fam_df.columns=["ukbb_FID", "ukbb_IID","ukbb_PID","ukbb_MID","ukbb_sex","ukbb_disease"]

common_samples=pd.merge(ukbb500k_fam_df,UKB_PPP_fam_df,left_on=["ukbb_FID", "ukbb_IID","ukbb_PID","ukbb_MID","ukbb_sex",
                            "ukbb_disease"], right_on=["ukpp_FID", "ukpp_IID","ukpp_PID","ukpp_MID","ukpp_sex","ukpp_disease"]).iloc[:,:2]
common_samples.to_csv("UKB_PPP_Samples_to_remove_removefrom_CogGWAS.tsv",index=None,header=None)

ukbb500k_UKB_PPP_df=pd.merge(ukbb500k_fam_df,UKB_PPP_fam_df,left_on=["ukbb_FID", "ukbb_IID","ukbb_PID","ukbb_MID","ukbb_sex",
                            "ukbb_disease"], right_on=["ukpp_FID", "ukpp_IID","ukpp_PID","ukpp_MID","ukpp_sex","ukpp_disease"],how="outer")
ukbb500k_UKB_specific_df=ukbb500k_UKB_PPP_df[ukbb500k_UKB_PPP_df["ukpp_IID"].isna()].iloc[:,0:2]


##Identify individula with cognitive phenotypes that are not part of 
cog_phenotype_df=pd.read_csv(cog_phenotype,sep=" ")
cog_phenotype_withoutUKB_PPP_df=pd.merge(ukbb500k_UKB_specific_df,cog_phenotype_df,left_on=["ukbb_FID","ukbb_IID"],
                                          right_on=["FID","IID"]).drop(["ukbb_FID","ukbb_IID"],axis=1)


cog_phenotype_withoutUKB_PPP_df.to_csv("ukbb500k_withoutUKB_PPP.fiall_baseline_model2.EUR.regenie.v3.txt",sep=" ",index=None)


# QC UKBB genotype data for Regenie based on the following code
os.system(f'''plink2 --threads 48 --pfile {geenotype_file} \
            --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --remove UKB_PPP_Samples_to_remove_removefrom_CogGWAS.tsv \
            --keep /edgehpc/dept/compbio/human_genetics/users/cchen11/ukbb_cogexome/program_v5/cvas_assoc_run2/ukbb500k.cogphenopc.EUR.feid.v1.keep \
            --write-snplist --write-samples --no-id-header --out ukb_chrauto_v2_snpqc.regenie_qc_pass.EUR ''')
            

#-------------------
# Regenie Step 1: fitting the null genome-wide ridge regression
# model 2 : sex + age + PC + assessment center
#-------------------

os.system(f'''regenie \
        --step 1 \
        --pgen /edgehpc/dept/compbio/human_genetics/ukbb_temp/Release2/Genotyped/ukb_chrauto_v2_snpqc \
        --extract ukb_chrauto_v2_snpqc.regenie_qc_pass.EUR.snplist \
        --keep ukb_chrauto_v2_snpqc.regenie_qc_pass.EUR.id \
        --covarFile ukbb500k_withoutUKB_PPP.fiall_baseline_model2.EUR.regenie.v3.txt \
        --phenoFile ukbb500k_withoutUKB_PPP.fiall_baseline_model2.EUR.regenie.v3.txt \
        --phenoColList fiall_baseline_irnt \
        --covarColList sex,fiall_baseline_age_mean0,fiall_baseline_age_mean02,fiall_baseline_age_mean0_sex,fiall_baseline_age_mean02_sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,assess_center10003,assess_center11001,assess_center11002,assess_center11003,assess_center11004,assess_center11005,assess_center11006,assess_center11007,assess_center11008,assess_center11009,assess_center11011,assess_center11012,assess_center11013,assess_center11014,assess_center11016,assess_center11017,assess_center11018,assess_center11020,assess_center11021,assess_center11022,assess_center11023 \
        --bsize 1000 \
        --threads 48 \
        --cv 10 \
        --out regenie_ukb_step1.fiall_baseline_irnt_model2''')

#-------------------
# Regenie Step 2: performing single-variant association tests
#-------------------

for chr in range(1,23):
   os.system(f''' 
echo """#!/bin/bash
#SBATCH --ntasks=4
#SBATCH --mem=260GB
#SBATCH --partition=cpu
#SBATCH --time=200:00:00
#SBATCH --job-name=regenie_ukb_step2_linear_model2_chr{chr}
#SBATCH --output=regenie_ukb_step2_linear_model2_chr{chr}.log
#SBATCH --error=regenie_ukb_step2_linear_model2_chr{chr}_err.log

    module load PLINK/2.00a2.3_x86_64
    module load Regenie/3.2.4-GCC-11.2.0
    module load BGEN-enkre/1.1.7-GCC-11.2.0
               
    regenie \
    --step 2 \
    --bgen /edgehpc/dept/compbio/human_genetics/ukbb_temp/Release2/Imputed/ukb_imp_chr{chr}_v3.bgen \
    --sample /edgehpc/dept/compbio/human_genetics/ukbb_temp/Release2/Imputed/ukb_imp_v3.sample \
    --keep ukb_chrauto_v2_snpqc.regenie_qc_pass.EUR.id \
    --covarFile ukbb500k_withoutUKB_PPP.fiall_baseline_model2.EUR.regenie.v3.txt \
    --phenoFile ukbb500k_withoutUKB_PPP.fiall_baseline_model2.EUR.regenie.v3.txt \
    --phenoColList fiall_baseline_irnt \
    --covarColList sex,fiall_baseline_age_mean0,fiall_baseline_age_mean02,fiall_baseline_age_mean0_sex,fiall_baseline_age_mean02_sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,assess_center10003,assess_center11001,assess_center11002,assess_center11003,assess_center11004,assess_center11005,assess_center11006,assess_center11007,assess_center11008,assess_center11009,assess_center11011,assess_center11012,assess_center11013,assess_center11014,assess_center11016,assess_center11017,assess_center11018,assess_center11020,assess_center11021,assess_center11022,assess_center11023 \
    --bsize 1000  \
    --pred regenie_ukb_step1.fiall_baseline_irnt_model2_pred.list \
    --chr {chr} --threads 48 \
    --out regenie_ukb_step2_linear_model2_chr{chr} """ >regenie_ukb_step2_linear_model2_chr{chr}.sh
    
    sbatch regenie_ukb_step2_linear_model2_chr{chr}.sh ''')
 

------------
# collect single-variant association test results
#-------------------

os.system(f''' cp regenie_ukb_step2_linear_model2_chr1_fiall_baseline_irnt.regenie regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.txt''')

for chr in range(2,23):
    os.system(f"""awk 'FNR>1 {{print $0}}' regenie_ukb_step2_linear_model2_chr{chr}_fiall_baseline_irnt.regenie >> regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.txt""")


#df=pd.read_csv("regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.txt",sep=" ")

os.system("gzip regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.txt")
os.system("cp regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.txt.gz /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/Cognition/Cognition_UKB/")




###### Filtering will be done during sumstat to vcf conversion time

#os.system(f"""echo "CHR BP SNP A1 A2 FREQA1 INFO BETA SE LOG10P P" > regenie_ukb_step2_linear_model2_fiall_baseline.regenie.info6maf005.summary.txt""")
#os.system(f"""zcat regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary.txt.gz | \
#              awk 'FNR>1{{if($6>=0.005 && $6<0.995 && $7>0.6) print $1,$2,$3,$5,$4,$6,$7,$9,$10,$12,10^(-1*$12)}}' >> \
#              regenie_ukb_step2_linear_model2_fiall_baseline.regenie.info6maf005.summary.txt""")


################ PC calculation
#plink2 \
#    --bfile /edgehpc/dept/compbio/human_genetics/users/cchen11/ukbb_gwas/pop_assign_v2/ukbb_gwas/pop_assign_v2/ukb_imp_chrauto_info08maf001_downsample1KGpruned \
#    --keep /edgehpc/dept/compbio/human_genetics/users/cchen11/ukbb_gwas/pop_assign_v2/ukbb_gwas/pop_assign_v2/ukb_imp_superpop_EUR.tsv \
#    --remove /edgehpc/dept/compbio/human_genetics/users/cchen11/ukbb_gwas/pop_assign_v2/ukbb_gwas/pop_assign_v2/ukb_imp_EUR_run1.outlier.remove \
#        /edgehpc/dept/compbio/human_genetics/users/cchen11/ukbb_gwas/pop_assign_v2/ukbb_gwas/pop_assign_v2/ukb_imp_EUR_run2.outlier.remove \
#            /edgehpc/dept/compbio/human_genetics/users/cchen11/ukbb_gwas/pop_assign_v2/ukbb_gwas/pop_assign_v2/ukb_imp_EUR_run3.outlier.remove \
#                --pca 40 approx --out --threads 24


