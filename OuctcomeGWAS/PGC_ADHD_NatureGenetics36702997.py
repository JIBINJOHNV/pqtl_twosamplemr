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

##File location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/PGC_ADHD
os.system("wget https://figshare.com/ndownloader/files/40036684")
os.system("mv 40036684 ADHD2022_iPSYCH_deCODE_PGC.meta.gz")
os.system('zgrep -v "##" ADHD2022_iPSYCH_deCODE_PGC.meta.gz  > ADHD2022_iPSYCH_deCODE_PGC.meta.tsv')

os.system("wget https://figshare.com/ndownloader/files/40036699")
os.system("mv 40036699 ADHD2022_README.pdf")



filename="ADHD2022_iPSYCH_deCODE_PGC.meta.tsv" #location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/PGC3_SCZ/UniqID
fdf=pd.read_csv(filename,sep=" ")

fdf3=fdf[['CHR', 'SNP', 'BP', 'A1', 'A2','FRQ_U_186843', 'INFO', 'OR', 'SE', 'P', 'Nca', 'Nco']]
fdf3[['BP','Nca', 'Nco']]=fdf3[['BP','Nca', 'Nco']].astype("int")
fdf3[['FRQ_U_186843', 'INFO', 'OR', 'SE', 'P']]=fdf3[['FRQ_U_186843', 'INFO', 'OR', 'SE', 'P']].astype("float")
fdf3=fdf3[['CHR','BP','SNP','A1','A2',"OR",'SE','Nco','Nca','P','FRQ_U_186843','INFO' ]]
fdf3.to_csv(f'{filename[:-4]}.tsv',sep="\t",index=None)


paramsdict={"chr_col": 0,
    "pos_col": 1,
    "snp_col": 2,
    "ea_col": 3,
    "oa_col": 4,
    "beta_col": 5,
    "se_col": 6,
    "ncontrol_col": 7,
    "ncase_col":8,
    "pval_col": 9,
    "eaf_col": 10,
    "imp_info_col":11,
    "delimiter": "\t",
    "header": "true",
    "build": "GRCh37"
    }

with open('ADHD2022_iPSYCH_deCODE_PGC_params.txt', 'w') as f:
  json.dump(paramsdict, f)


##Convert to vcf format 
path=os.getcwd()
ID="PGC_ADHD2022_iPSYCH_deCODE"

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename[:-4]}.vcf \
    --data {path}/{filename[:-4]}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json ADHD2022_iPSYCH_deCODE_PGC_params.txt \
    --id {ID} > {filename[:-4]}.error 2>&1 ''' )


##Convert Grch37 to GRCh38
os.system(f'''CrossMap.py vcf /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GRCh37_to_GRCh38.chain {path}/{filename[:-4]}.vcf.gz \
            /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa {path}/{filename[:-4]}_GRCh38.vcf''')

##Sort vcf file
os.system(f"bcftools sort {path}/{filename[:-4]}_GRCh38.vcf | bgzip -c > {path}/{filename[:-4]}_GRCh38.vcf.gz")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38.vcf.gz')
os.system(f"rm {path}/{filename[:-4]}_GRCh38.vcf ")


##Ad dbsnp 156 to the file
os.system(f'''bcftools annotate --threads  10 \
    -a /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GCF_000001405_cleaed.40_Msplited.vcf.gz -c ID \
    -o {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz -O z {path}/{filename[:-4]}_GRCh38.vcf.gz''')

os.system(f'tabix -f -p vcf {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz')
os.system(f"rm {path}/{filename[:-4]}_GRCh38.vcf.gz* {path}/{filename[:-4]}.vcf.gz")

##Ading dbsnp 156 to the format section
os.system(f'zgrep -v "##" {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz > {path}/{filename[:-4]}_GRCh38_rsid156.tab')
os.system(f'zgrep  "##" {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz > {path}/{filename[:-4]}_GRCh38_rsid156_header.txt')




df=pd.read_csv(f"{path}/{filename[:-4]}_GRCh38_rsid156.tab",sep="\t")
format_df2=df["FORMAT"].str.split(":",expand=True)
format_df2[7]="ID"
df["FORMAT"]=format_df2[0]+":"+format_df2[1]+":"+format_df2[2]+":"+format_df2[3]+":"+format_df2[4]+":"+format_df2[5]+":"+format_df2[6]+":"+format_df2[7]

format_df=df[ID].str.split(":",expand=True)
df["ID"]=np.where(df["ID"].str.contains("rs"),df["ID"],df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"])
format_df[7]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]+":"+format_df[6]+":"+format_df[7]

df["FORMAT"].str.split(":",expand=True).isna().sum()
df[ID].str.split(":",expand=True).isna().sum()

df.to_csv(f"{path}/{filename[:-4]}_GRCh38_rsid156.tab",sep="\t",index=None)

os.system(f"cat {path}/{filename[:-4]}_GRCh38_rsid156_header.txt {path}/{filename[:-4]}_GRCh38_rsid156.tab | bgzip -c > {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz')



##Ading UniqID to the format section
df["ID"]=df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"]
format_df[7]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]+":"+format_df[6]+":"+format_df[7]

df.to_csv(f"{path}/{filename[:-4]}_GRCh38_UniqID.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-4]}_GRCh38_rsid156_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab | bgzip -c > {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz')

os.system(f" rm {path}/{filename[:-4]}_GRCh38_rsid156_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab {path}/{filename[:-4]}_GRCh38_rsid156.tab  ")



          





