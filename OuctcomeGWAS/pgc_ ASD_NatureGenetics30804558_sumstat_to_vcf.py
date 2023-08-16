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
fdf3["Beta"]=fdf3["OR"]
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














