
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
filename=args.genomeDir

#filename="Sample_list_Part1.tsv"

fdf=pd.read_csv(filename,sep="\t")

path=os.getcwd()

n=1
for file in fdf['SampleNames'].unique():
    print(f'''Now running {n} {file.split("/")[-1][:-3]}''')
    df=pd.read_csv(file,sep="\t")
    df=df[df['LOG10P']>=5]
    df['LOG10P']=np.where(df['LOG10P']<200,df['LOG10P'],199)
    df=df[['CHROM', 'GENPOS','ID','ALLELE1','ALLELE0','BETA','SE','N' ,'LOG10P','A1FREQ','INFO','GeneSymbol']]
    df['LOG10P']=np.power(10,-df['LOG10P'])
    ID=df.iloc[0,11]
    df['ID']=df['ID'].str.replace(":","_")
    df['CHROM']=np.where(df['CHROM']==23,"X",df['CHROM'])
    df.to_csv(f'{file.split("/")[-1][:-3]}',sep="\t",index=None)
    n=n+1
    os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Discovery_Genewise_Merged/gwas2vcf/main.py \
        --out {path}/{file.split("/")[-1][:-7]}.vcf \
        --data {path}/{file.split("/")[-1][:-3]} \
        --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        --dbsnp /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GCF_000001405_cleaed.40_Msplited.vcf.gz \
        --json params.txt \
        --id {ID} ''' )
    os.system(f'tabix -f -p vcf  {path}/{file.split("/")[-1][:-7]}.vcf.gz')

#paramsdict={"chr_col": 0,
#    "pos_col": 1,
#    "snp_col": 2,
#    "ea_col": 3,
#    "oa_col": 4,
#    "beta_col": 5,
#    "se_col": 6,
#    "ncontrol_col": 7,
#    "pval_col": 8,
#    "eaf_col": 9,
#    "imp_info_col":10,
#    "delimiter": "\t",
#    "header": "true",
#    "build": "GRCh38"
#    }
#
#with open('params.txt', 'w') as f:
#  json.dump(paramsdict, f)

#paramsdict={"chr_col": 0,
#    "pos_col": 2,
#    "snp_col": 1,
#    "ea_col": 3,
#    "oa_col": 4,
#    "beta_col": 8,
#    "se_col": 9,
#    "ncontrol_col": 12,
#    "ncase_col":11,
#    "pval_col": 10,
#    "eaf_col": 5,
#    "imp_info_col":7,
#    "delimiter": "\t",
#    "header": "true",
#    "build": "GRCh37"}



CrossMap.py vcf --compress /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GRCh37_to_GRCh38.chain PGC3_SCZ_wave3.european.autosome.public.v3.vcf.gz \
            /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa out_vcf


bcftools annotate --threads  10 \
    -a /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GCF_000001405_cleaed.40_Msplited.vcf.gz -c ID \
   -o PGC3_SCZ_wave3.european.autosome.public.v3_hg38_rsid156.vcf.gz \
   -O z PGC3_SCZ_wave3.european.autosome.public.v3_hg38.vcf.gz 


