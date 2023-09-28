

##Before running run below mentioned unix command
#Conda deactivate
module load Python/3.8.6-GCCcore-10.2.0
module load libxml2/2.9.8-GCCcore-6.4.0
module load OpenSSL/1.1
module load R/4.1.3-foss-2021b
module load GMP/6.2.1-GCCcore-11.2.0
module load BCFtools/1.14-GCC-11.2.0

#source /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/env/bin/activate


import pandas as pd
import numpy as np
import tempfile
import json,os
import argparse


pd.set_option('display.float_format', '{:.2E}'.format)
print(tempfile.gettempdir())

##File location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/Cognition_Savage_iq
## The file shared by MAx
os.system('zgrep -v "##" Cognition_Savage_iq_metaanalysis_noUKB.out.gz > Cognition_Savage_iq_metaanalysis_noUKB.tsv')


#parser=argparse.ArgumentParser(description="It is for Incorporating additional information and Initial filtering of Annovar out put file")
#parser.add_argument('-filename','--filename', help="file name ", required=True)
#args=parser.parse_args()
#filename=args.filename

filename="Cognition_Savage_iq_metaanalysis_noUKB.tsv" #location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/Cognition_Savage_iq/"
fdf=pd.read_csv(filename,sep="\t")
fdf=fdf[['CHR','SNP','BP', 'Allele1', 'Allele2','totalN','Zscore', 'P', 'EAF_HRC']]

#https://ctg.cncr.nl/documents/p1651/readme.txt ;  https://www.biostars.org/p/319584/
fdf['beta'] = fdf['Zscore'] / np.sqrt(2 * fdf['EAF_HRC'] * (1 - fdf['EAF_HRC']) * (fdf['totalN'] + fdf['Zscore']**2))
fdf['se']   = 1 / np.sqrt(2 * fdf['EAF_HRC'] * (1 - fdf['EAF_HRC']) * (fdf['totalN'] + fdf['Zscore']**2))


fdf=fdf[['CHR','SNP','BP', 'Allele1','Allele2','EAF_HRC','beta','se','P','totalN']]
fdf['Allele1']=fdf['Allele1'].str.upper()
fdf['Allele2']=fdf['Allele2'].str.upper()

fdf3=fdf[~fdf["CHR"].isna()]
fdf3["CHR"]=fdf3["CHR"].astype("int").astype("str")
fdf3[['BP','totalN']]=fdf3[['BP','totalN']].astype("int")

fdf3[['EAF_HRC', 'beta','se','P']]=fdf3[['EAF_HRC', 'beta','se','P']].astype("float")
fdf3=fdf3[['CHR','BP','SNP','Allele1','Allele2','beta','se','totalN','P','EAF_HRC']]
fdf3.to_csv(f'{filename[:-4]}.tsv',sep="\t",index=None)

paramsdict={"chr_col": 0,
    "pos_col": 1,
    "snp_col": 2,
    "ea_col": 3,
    "oa_col": 4,
    "beta_col": 5,
    "se_col": 6,
    "ncontrol_col": 7,
    "pval_col": 8,
    "eaf_col": 9,
    "delimiter": "\t",
    "header": "true",
    "build": "GRCh37"
    }

with open('Cognition_Savage_iq_metaanalysis_noUKB_params.txt', 'w') as f:
  json.dump(paramsdict, f)

##Convert to vcf format 
path=os.getcwd()
ID="Cognition_Savage_NoUKB"

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename[:-4]}.vcf \
    --data {path}/{filename[:-4]}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json Cognition_Savage_iq_metaanalysis_noUKB_params.txt \
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
format_df2[5]="ID"
df["FORMAT"]=format_df2[0]+":"+format_df2[1]+":"+format_df2[2]+":"+format_df2[3]+":"+format_df2[4]+":"+format_df2[5]

format_df=df[ID].str.split(":",expand=True)
df["ID"]=np.where(df["ID"].str.contains("rs"),df["ID"],df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"])
format_df[5]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[4]

df["FORMAT"].str.split(":",expand=True).isna().sum()
df[ID].str.split(":",expand=True).isna().sum()

df.to_csv(f"{path}/{filename[:-4]}_GRCh38_rsid156.tab",sep="\t",index=None)

os.system(f"cat {path}/{filename[:-4]}_GRCh38_rsid156_header.txt {path}/{filename[:-4]}_GRCh38_rsid156.tab | bgzip -c > {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz')



##Ading UniqID to the format section
df["ID"]=df["#CHROM"].astype(str)+"_"+df["POS"].astype(str)+"_"+df["REF"]+"_"+df["ALT"]
format_df[5]=df["ID"]
format_df=format_df.astype("str")
df[ID]=format_df[0]+":"+format_df[1]+":"+format_df[2]+":"+format_df[3]+":"+format_df[4]+":"+format_df[5]

df.to_csv(f"{path}/{filename[:-4]}_GRCh38_UniqID.tab",sep="\t",index=None)
os.system(f"cat {path}/{filename[:-4]}_GRCh38_rsid156_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab | bgzip -c > {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz ")
os.system(f'tabix -f -p vcf  {path}/{filename[:-4]}_GRCh38_UniqID.vcf.gz')

df2=pd.concat([df.iloc[:,0:5],format_df],axis=1).iloc[:,:-1]
df2.to_csv(f"{path}/{filename[:-4]}_GRCh38_UniqID_ForMeta.tsv",sep="\t",index=None)


os.system(f" rm {path}/{filename[:-4]}_GRCh38_rsid156_header.txt {path}/{filename[:-4]}_GRCh38_UniqID.tab {path}/{filename[:-4]}_GRCh38_rsid156.tab  ")

