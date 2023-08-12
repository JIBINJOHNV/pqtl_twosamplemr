##Before running run below mentioned unix command
#Conda deactivate
#module load Python/3.8.6-GCCcore-10.2.0
#module load libxml2/2.9.8-GCCcore-6.4.0
#module load OpenSSL/1.1
#module load R/4.1.3-foss-2021b
#module load GMP/6.2.1-GCCcore-11.2.0
#source /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/env/bin/activate


import pandas as pd
import numpy as np
import tempfile
import json,os
import argparse


pd.set_option('display.float_format', '{:.2E}'.format)
print(tempfile.gettempdir())


parser=argparse.ArgumentParser(description="It is for Incorporating additional information and Initial filtering of Annovar out put file")
parser.add_argument('-filename','--filename', help="file name ", required=True)


args=parser.parse_args()
filename=args.filename

filename="daner_bip_pgc3_nm_noukbiobank.gz" #location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/bip_pgc3

fdf=pd.read_csv(filename,sep="\t")
fdf2=fdf[['CHR', 'SNP', 'BP', 'A1', 'A2','FRQ_U_313436', 'INFO','OR', 'SE', 'P','Nca', 'Nco']]
fdf3=fdf2[~fdf2["SNP"].isna()]
fdf3[['BP','Nca', 'Nco']]=fdf3[['BP','Nca', 'Nco']].astype("int")
fdf3[['FRQ_U_313436', 'INFO', 'OR', 'SE', 'P']]=fdf3[['FRQ_U_313436', 'INFO', 'OR', 'SE', 'P']].astype("float")
fdf3["Beta"]=np.log(fdf3["OR"])
fdf3.drop("OR",axis=1,inplace=True)
fdf3=fdf3[['CHR','BP','SNP','A1','A2',"Beta",'SE','Nco','Nca','P','FRQ_U_313436','INFO' ]]
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
    "eaf_col": 10,
    "imp_info_col":11,
    "delimiter": "\t",
    "header": "true",
    "build": "GRCh38"
    }

with open('daner_bip_pgc3_params.txt', 'w') as f:
  json.dump(paramsdict, f)



path=os.getcwd()

ID="BIP_PGC3_noukb"

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename}[:-3].vcf \
    --data {path}/{filename[:-3]}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json daner_bip_pgc3_params.txt \
    --id {ID} ''' )

os.system(f'tabix -f -p vcf  {path}/{filename}[:-3].vcf.gz')
