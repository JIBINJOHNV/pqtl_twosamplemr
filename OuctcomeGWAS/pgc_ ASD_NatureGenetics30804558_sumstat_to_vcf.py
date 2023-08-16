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

##File location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/PGC_ASD
os.system("wget https://figshare.com/ndownloader/files/28169292")
os.system("mv 28169292 iPSYCH-PGC_ASD_Nov2017.gz")

os.system("wget https://figshare.com/ndownloader/files/28169289")
os.system("mv 28169289 iPSYCH-PGC_ASD_Nov2017_readme.pdf")




filename="iPSYCH-PGC_ASD_Nov2017.gz" #location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/bip_pgc3
fdf=pd.read_csv(filename,sep="\t")
fdf['Nca']=18381
fdf['Nco']=27969

fdf2=fdf[['CHR', 'SNP', 'BP', 'A1', 'A2','INFO','OR', 'SE', 'P','Nca', 'Nco']]
fdf3=fdf2[~fdf2["SNP"].isna()]
fdf3[['BP','Nca', 'Nco']]=fdf3[['BP','Nca', 'Nco']].astype("int")
fdf3[['INFO', 'OR', 'SE', 'P']]=fdf3[['INFO', 'OR', 'SE', 'P']].astype("float")
fdf3["Beta"]=np.log(fdf3["OR"])
fdf3.sort_values(by="Beta")
fdf3.drop("OR",axis=1,inplace=True)
fdf3=fdf3[['CHR','BP','SNP','A1','A2',"Beta",'SE','Nco','Nca','P','INFO' ]]
fdf3.to_csv(f'{filename[:-3]}.tsv',sep="\t",index=None)


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
    "imp_info_col":10,
    "delimiter": "\t",
    "header": "true",
    "build": "GRCh37"
    }

with open('asd_pgc3_params.txt', 'w') as f:
  json.dump(paramsdict, f)


path=os.getcwd()
ID="ASD_PGC"

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename[:-3]}.vcf \
    --data {path}/{filename[:-3]}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json asd_pgc3_params.txt \
    --id {ID} > {filename[:-3]}.error 2>&1 ''' )

#os.system(f'tabix -f -p vcf  {path}/{filename[:-3]}.vcf.gz')


os.system(f'''CrossMap.py vcf /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GRCh37_to_GRCh38.chain {path}/{filename[:-3]}.vcf.gz \
            /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa {path}/{filename[:-3]}_GRCh38.vcf''')

os.system(f"bcftools sort {path}/{filename[:-3]}_GRCh38.vcf | bgzip -c > {path}/{filename[:-3]}_GRCh38.vcf.gz")
os.system(f'tabix -f -p vcf  {path}/{filename[:-3]}_GRCh38.vcf.gz')
os.system("rm daner_bip_pgc3_nm_noukbiobank_GRCh38.vcf ")

os.system(f'''bcftools annotate --threads  10 \
    -a /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GCF_000001405_cleaed.40_Msplited.vcf.gz -c ID \
    -o {path}/{filename[:-3]}_GRCh38_rsid156.vcf.gz -O z {path}/{filename[:-3]}_GRCh38.vcf.gz''')

os.system(f'tabix -f -p vcf  {path}/{filename[:-3]}_GRCh38_rsid156.vcf.gz')
os.system(f'zgrep -v "##" {path}/{filename[:-3]}_GRCh38_rsid156.vcf.gz > {path}/{filename[:-3]}_GRCh38_rsid156.tab')
os.system(f'zgrep  "##" {path}/{filename[:-3]}_GRCh38_rsid156.vcf.gz > {path}/{filename[:-3]}_GRCh38_rsid156_header.txt')


df=pd.read_csv(f"{path}/{filename[:-3]}_GRCh38_rsid156.tab",sep="\t")
format_df=df[ID].str.split(":",expand=True)
df["ID"]=np.where(df["ID"].str.contains("rs"),df["ID"],df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"])
format_df[6]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]+":"+format_df[6]

format_df2=df["FORMAT"].str.split(":",expand=True)
format_df2[6]="ID"
df['FORMAT']=format_df2[0]+":"+format_df2[1]+":"+format_df2[2]+":"+format_df2[3]+":"+format_df2[4]+":"+format_df2[5]+":"+format_df2[6]

df['FORMAT'].str.split(":",expand=True).isna().sum()
df[ID].str.split(":",expand=True).isna().sum()
print(df)

df.to_csv(f"{path}/{filename[:-3]}_GRCh38_rsid156.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-3]}_GRCh38_rsid156_header.txt {path}/{filename[:-3]}_GRCh38_rsid156.tab | bgzip -c > {path}/{filename[:-3]}_GRCh38_rsid156.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-3]}_GRCh38_rsid156.vcf.gz')


##Replace rsid with Uniqid
df["ID"]=df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"]
format_df[6]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]+":"+format_df[6]

df['FORMAT'].str.split(":",expand=True).isna().sum()
df[ID].str.split(":",expand=True).isna().sum()
print(df)

df.to_csv(f"{path}/{filename[:-3]}_GRCh38_UniqID.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-3]}_GRCh38_rsid156_header.txt {path}/{filename[:-3]}_GRCh38_UniqID.tab | bgzip -c > {path}/{filename[:-3]}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-3]}_GRCh38_UniqID.vcf.gz')
os.system(f" rm {path}/{filename[:-3]}_GRCh38_rsid156_header.txt {path}/{filename[:-3]}_GRCh38_UniqID.tab {path}/{filename[:-3]}_GRCh38_rsid156.tab  ")
os.system(f"rm {path}/{filename[:-3]}.vcf.gz* {path}/{filename[:-3]}_GRCh38.vcf.gz*")

