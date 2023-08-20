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



os.system('zgrep -v "##" pgcAN2.2019-07.vcf.tsv.gz  > pgcAN2.2019-07.vcf.tsv')
filename="pgcAN2.2019-07.vcf.tsv" #location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/PGC_EatingDisorder
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


os.system(f'''bcftools annotate --threads  10 \
    -a /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GCF_000001405_cleaed.40_Msplited.vcf.gz -c ID \
    -o {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz -O z {path}/{filename[:-4]}_GRCh38.vcf.gz''')


os.system(f'''bcftools annotate --threads  10 --collapse none \
              -a /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/1000GENOMES-phase_3_GRCh38_Msplited.vcf.gz -c 'INFO/EUR:=ID' \
              {path}/{filename[:-4]}_GRCh38_rsid156_1kgEUR_VAF.vcf.gz  -O z {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz''')



##Since allele frequency not present we useed 1000 gENOME EUR Freq ; location /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/
os.system("wget https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz")
os.system("""bcftools norm -m-any --check-ref -w \
              > -f /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
              > /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/1000GENOMES-phase_3.vcf.gz | bgzip -c > \
              /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/1000GENOMES-phase_3_GRCh38_Msplited.vcf.gz""")

os.system("tabix -p vcf 1000GENOMES-phase_3_GRCh38_Msplited.vcf.gz")










