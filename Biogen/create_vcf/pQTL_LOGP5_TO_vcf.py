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
