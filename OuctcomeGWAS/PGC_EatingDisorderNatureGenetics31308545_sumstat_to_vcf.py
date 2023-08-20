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
fdf3[['FRQ_U_186843', 'INFO', 'OR', 'SE', 'P']]=fdf3[['FRQ_U_186843', 'INFO', 'OR', 'SE', 'P']].astype("float")
fdf3=fdf3[['CHR','BP','SNP','A1','A2',"OR",'SE','Nco','Nca','P','FRQ_U_186843','INFO' ]]
fdf3[fdf3["OR"]<0]





