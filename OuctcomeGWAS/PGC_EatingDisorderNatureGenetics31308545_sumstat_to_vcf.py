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

##File location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/PGC_EatingDisorder
os.system("wget https://figshare.com/ndownloader/files/28169271")
os.system("mv 28169271 pgcAN2.2019-07.vcf.tsv.gz")
os.system("wget https://figshare.com/ndownloader/files/28169268")
os.system("mv 28169268 an2019.readme.pdf")



os.system('zgrep -v "##" pgcAN2.2019-07.vcf.tsv.gz  > pgcAN2.2019-07.tsv')
filename="pgcAN2.2019-07.tsv" #location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/PGC_EatingDisorder
fdf=pd.read_csv(filename,sep="\t")


fdf3=fdf[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'BETA', 'SE', 'PVAL','IMPINFO','NCAS', 'NCON']]
fdf3[['POS','NCAS', 'NCON']]=fdf3[['POS','NCAS', 'NCON']].astype("int")
fdf3[['IMPINFO', 'BETA', 'SE', 'PVAL']]=fdf3[['IMPINFO', 'BETA', 'SE', 'PVAL']].astype("float")
fdf3=fdf3[['CHROM','POS','ID','ALT','REF',"BETA",'SE','NCON','NCAS','PVAL','IMPINFO' ]]
fdf3[fdf3["BETA"]<0]


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
    "imp_info_col":10,
    "delimiter": "\t",
    "header": "true",
    "build": "GRCh37"
    }

with open('pgcAN2.2019-07_params.txt', 'w') as f:
  json.dump(paramsdict, f)


##Convert to vcf format 
path=os.getcwd()
ID="PGC_AN2"

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename[:-4]}.vcf \
    --data {path}/{filename[:-4]}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json pgcAN2.2019-07_params.txt \
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


##Since allele frequency not present we useed 1000 gENOME EUR Freq ; location /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/
os.system("wget https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz")
os.system("""bcftools norm -m-any --check-ref -w \
              > -f /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
              > /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/1000GENOMES-phase_3.vcf.gz | bgzip -c > \
              /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/1000GENOMES-phase_3_GRCh38_Msplited.vcf.gz""")

os.system("tabix -p vcf 1000GENOMES-phase_3_GRCh38_Msplited.vcf.gz")


##Adding Europe 1kg allele frequency
os.system(f'''bcftools annotate --threads  10 --collapse none \
              -a /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/1000GENOMES-phase_3_GRCh38_Msplited.vcf.gz -c 'INFO/EUR' \
              -o {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.vcf.gz  -O z {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz''')

os.system(f"rm  {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz*  {path}/{filename[:-4]}_GRCh38.vcf.gz* {path}/{filename[:-4]}.vcf*")



##Ading dbsnp 156 to the format section
os.system(f'zgrep -v "##" {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.vcf.gz > {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.tab')
os.system(f'zgrep  "##" {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.vcf.gz > {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF_header.txt')


df=pd.read_csv(f"{path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.tab",sep="\t")
df=df[df["INFO"].str.contains("EUR")]

format_df=df[ID].str.split(":",expand=True)
df["ID"]=np.where(df["ID"].str.contains("rs"),df["ID"],df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"])
df["INFO"]=df["INFO"].str.replace("EUR","AF")
format_df[6]=df["ID"]
format_df[7]=df["INFO"].str.replace("AF=","")
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]+":"+format_df[6]+":"+format_df[7]

format_df2=df["FORMAT"].str.split(":",expand=True)
format_df2[6]="ID"
format_df2[7]="AF"
df['FORMAT']=format_df2[0]+":"+format_df2[1]+":"+format_df2[2]+":"+format_df2[3]+":"+format_df2[4]+":"+format_df2[5]+":"+format_df2[6]+":"+format_df2[7]
df['FORMAT'].str.split(":",expand=True).isna().sum()
df[ID].str.split(":",expand=True).isna().sum()
print(df)

df.to_csv(f"{path}/{filename[:-3]}_GRCh38_rsid156_1kgEUR_VAF.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF_header.txt {path}/{filename[:-3]}_GRCh38_rsid156_1kgEUR_VAF.tab | bgzip -c > {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.gz')


##Replace rsid with Uniqid
df["ID"]=df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"]
format_df[6]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]+":"+format_df[6]+":"+format_df[7]

df['FORMAT'].str.split(":",expand=True).isna().sum()
df[ID].str.split(":",expand=True).isna().sum()
print(df)

df.to_csv(f"{path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF_header.txt {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.tab | bgzip -c > {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz')
os.system(f" rm {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF_header.txt {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.tab {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.tab  ")
os.system(f"rm {path}/{filename[:-4]}.vcf.gz* {path}/{filename[:-4]}_GRCh38.vcf.gz*")




