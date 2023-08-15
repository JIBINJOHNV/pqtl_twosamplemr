
##Before running run below mentioned unix command
#Conda deactivate
#module load Python/3.8.6-GCCcore-10.2.0
#module load libxml2/2.9.8-GCCcore-6.4.0
#module load OpenSSL/1.1
#module load R/4.1.3-foss-2021b
#module load GMP/6.2.1-GCCcore-11.2.0
module load BCFtools/1.14-GCC-11.2.0

#source /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/env/bin/activate


import pandas as pd
import numpy as np
import tempfile
import json,os
import argparse


pd.set_option('display.float_format', '{:.2E}'.format)
print(tempfile.gettempdir())

##File location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/MDD_PGC_Depression
os.system("wget https://datashare.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt")
os.system("wget https://datashare.ed.ac.uk/bitstream/handle/10283/3203/ReadMe.txt")
#The file PGC_UKB_depression_genome-wide.txt contains genome-wide summary statistics from a meta-analysis of the 33 cohorts of the Psychiatric Genomics Consortium 
#(excluding UK Biobank and 23andMe data) as described in Wray et al. (2018), https://doi.org/10.1038/s41588-018-0090-3 and the broad depression phenotype in the full release of the UK Biobank as described in Howard et al. (2018), https://doi.org/10.1038/s41467-018-03819-3.
#The total number of individuals in this data is 500,199 (170,756 cases and 329,443 controls) with 8,483,301 variants analysed.
#Saples	Cases	Controls
#PGC2913	16823	25632
#deCODE13	1980	9536
#GenScotland	997	6358
#GERA	7162	38307
#iPSYCH	18629	17841
#Totall	45591	97674

#parser=argparse.ArgumentParser(description="It is for Incorporating additional information and Initial filtering of Annovar out put file")
#parser.add_argument('-filename','--filename', help="file name ", required=True)
#args=parser.parse_args()
#filename=args.filename


filename="PGC_UKB_depression_genome-wide.txt" #location /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/MDD_PGC_Depression
fdf=pd.read_csv(filename,sep=" ")
fdf.columns=['ID','A1', 'A2',"Freq",'BETA','SE','PVAL']
fdf[["Freq",'BETA','SE','PVAL']]=fdf[["Freq",'BETA','SE','PVAL']].astype("float")
fdf["A1"]=fdf["A1"].str.upper()
fdf["A2"]=fdf["A2"].str.upper()


##Sincee the GWAS file dont have Chromosome position we will extract that frm dbsnp 156 hg38 file
#os.system(f'''zgrep -v "##" GCF_000001405_cleaed.40_Msplited.vcf.gz | egrep "SNV|POS" | awk -F"\t" '{print $1,$2,$3,$4,$5}' > GCF_000001405_cleaed.40_Msplited_SNP_ChrPosRefAlt.txt''')

snp=pd.read_csv("/edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GCF_000001405_cleaed.40_Msplited_SNP_ChrPosRefAlt.txt",sep=" ")
merged=pd.merge(fdf,snp,on="ID",how="left")

merged.to_csv("GCF_000001405_cleaed.40_Msplited_SNP_ChrPosRefAlt_PGC_UKB_depression_genome-wide.csv"index=None)
merged["#CHROM"]=merged["#CHROM"].astype("str")
merged2=merged[~merged["#CHROM"].isna()]
merged2=merged2[~merged2["#CHROM"].str.contains("N")]
merged2=merged2[merged2["#CHROM"]!='nan']

fdf_merged=merged2[ ((merged2["A1"]==merged2["REF"]) | (merged2["A1"]==merged2["ALT"])) & ((merged2["A2"]==merged2["REF"]) | (merged2["A2"]==merged2["ALT"]))  ]
fdf_merged.drop(["REF","ALT"],axis=1,inplace=True)
fdf_merged["NCAS"]=45591
fdf_merged["NCON"]=97674


fdf_merged=fdf_merged[['#CHROM','ID','POS','A1','A2','Freq','BETA', 'SE', 'PVAL','NCAS', 'NCON']]
fdf_merged.columns=['CHROM','ID','POS','A1','A2','FCON','BETA', 'SE', 'PVAL', 'NCAS', 'NCON']
fdf_merged=fdf_merged.sort_values(by=['CHROM','POS'])
fdf_merged[['CHROM','POS']]=fdf_merged[['CHROM','POS']].astype("int")
fdf_merged=fdf_merged[['CHROM','POS','ID','A1','A2','BETA','SE','NCON','NCAS','PVAL','FCON']]


##Convert to vcf format 
path=os.getcwd()
ID="PGC_MDD_Depression" 

fdf_merged.to_csv(f"{path}/{filename[:-4]}_GRCh38_ChrPosition.tab",sep="\t",index=None)

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
    "delimiter": "\t",
    "header": "true",
    "build": "GRCh38"
    }

with open('PGC_MDD_Depression_params.txt', 'w') as f:
  json.dump(paramsdict, f)


os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{filename[:-4]}.vcf \
    --data {path}/{filename[:-4]}_GRCh38_ChrPosition.tab \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --json PGC_MDD_Depression_params.txt \
    --id {ID} > {filename[:-4]}.error 2>&1 ''' )


##Ad dbsnp 156 to the file
os.system(f'''bcftools annotate --threads  10 \
    -a /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GCF_000001405_cleaed.40_Msplited.vcf.gz -c ID \
    -o {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz -O z {path}/{filename[:-4]}.vcf.gz''')

os.system(f'tabix -f -p vcf {path}/{filename[:-4]}_GRCh38_rsid156.vcf.gz')
os.system(f"rm  {path}/{filename[:-4]}.vcf.gz")

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








