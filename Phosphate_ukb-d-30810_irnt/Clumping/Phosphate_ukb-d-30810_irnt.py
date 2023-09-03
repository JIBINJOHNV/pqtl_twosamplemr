
##Before running run below mentioned unix command
#Conda deactivate
#module load Python/3.8.6-GCCcore-10.2.0
#module load libxml2/2.9.8-GCCcore-6.4.0
#module load OpenSSL/1.1
#module load R/4.1.3-foss-2021b
#module load GMP/6.2.1-GCCcore-11.2.0
#module load BCFtools/1.14-GCC-11.2.0

#source /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/env/bin/activate


import pandas as pd
import numpy as np
import tempfile
import json,os
import argparse


pd.set_option('display.float_format', '{:.2E}'.format)
print(tempfile.gettempdir())

##File location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/Phosphate_ukb-d-30810_irnt
os.system("wget https://gwas.mrcieu.ac.uk/files/ukb-d-30810_irnt/ukb-d-30810_irnt.vcf.gz")
os.system("wget https://gwas.mrcieu.ac.uk/files/ukb-d-30810_irnt/ukb-d-30810_irnt.vcf.gz.tbi")

filename="ukb-d-30810_irnt.vcf.gz"
path=os.getcwd()
ID="ukb-d-30810_irnt"
ID2="Phosphate_ukb-d-30810_irnt"

## The file is in vcf format and Grch37  so will convert to hg38



##Convert Grch37 to GRCh38
os.system(f'''CrossMap.py vcf /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GRCh37_to_GRCh38.chain {path}/{filename} \
            /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa {path}/{filename[:-7]}_GRCh38.vcf''')

##Sort vcf file
os.system(f"bcftools sort {path}/{filename[:-7]}_GRCh38.vcf | bgzip -c > {path}/{filename[:-7]}_GRCh38.vcf.gz")
os.system(f'tabix -f -p vcf  {path}/{filename[:-7]}_GRCh38.vcf.gz')
os.system(f"rm {path}/{filename[:-7]}_GRCh38.vcf ")


##Ad dbsnp 156 to the file
os.system(f'''bcftools annotate --threads  10 \
    -a /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GCF_000001405_cleaed.40_Msplited.vcf.gz -c ID \
    -o {path}/{filename[:-7]}_GRCh38_rsid156.vcf.gz -O z {path}/{filename[:-7]}_GRCh38.vcf.gz''')

os.system(f'tabix -f -p vcf {path}/{filename[:-7]}_GRCh38_rsid156.vcf.gz')
os.system(f"rm {path}/{filename[:-7]}_GRCh38.vcf.gz* {path}/{filename[:-7]}.vcf.gz")

##Ading dbsnp 156 to the format section
os.system(f'zgrep -v "##" {path}/{filename[:-7]}_GRCh38_rsid156.vcf.gz > {path}/{filename[:-7]}_GRCh38_rsid156.tab')
os.system(f'zgrep  "##" {path}/{filename[:-7]}_GRCh38_rsid156.vcf.gz > {path}/{filename[:-7]}_GRCh38_rsid156_header.txt')


df=pd.read_csv(f"{path}/{filename[:-7]}_GRCh38_rsid156.tab",sep="\t")
df=df[df["#CHROM"]!="X"]


format_df2=df["FORMAT"].str.split(":",expand=True)
format_df2[5]="ID"
df["FORMAT"]=format_df2[0]+":"+format_df2[1]+":"+format_df2[2]+":"+format_df2[3]+":"+format_df2[4]+":"+format_df2[5]

format_df=df[ID].str.split(":",expand=True)
df["ID"]=np.where(df["ID"].str.contains("rs"),df["ID"],df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"])
format_df[5]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]

df["FORMAT"].str.split(":",expand=True).isna().sum()
df[ID].str.split(":",expand=True).isna().sum()

df.rename(columns={ID:ID2},inplace=True)

df.to_csv(f"{path}/{filename[:-7]}_GRCh38_rsid156.tab",sep="\t",index=None)

os.system(f"cat {path}/{filename[:-7]}_GRCh38_rsid156_header.txt {path}/{filename[:-7]}_GRCh38_rsid156.tab | bgzip -c > {path}/{filename[:-7]}_GRCh38_rsid156.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-7]}_GRCh38_rsid156.vcf.gz')



##Ading UniqID to the format section
df["ID"]=df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"]
format_df[5]=df["ID"]
format_df=format_df.astype("str")
df[ID2]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]

df.to_csv(f"{path}/{filename[:-7]}_GRCh38_UniqID.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-7]}_GRCh38_rsid156_header.txt {path}/{filename[:-7]}_GRCh38_UniqID.tab | bgzip -c > {path}/{filename[:-7]}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-7]}_GRCh38_UniqID.vcf.gz')

os.system(f" rm {path}/{filename[:-7]}_GRCh38_rsid156_header.txt {path}/{filename[:-7]}_GRCh38_UniqID.tab {path}/{filename[:-7]}_GRCh38_rsid156.tab  ")

