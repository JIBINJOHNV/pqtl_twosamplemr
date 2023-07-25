#srun -n 12 --mem=124G --time=150:01:00  --pty bash
#module load Python/3.8.8-GCCcore-11.2.0
#module load libxml2/2.9.8-GCCcore-6.4.0
#module load OpenSSL/1.1
#source /edgehpc/dept/human_genetics/users/jjohn1/Discovery_Genewise_Merged/gwas2vcf/env/bin/activate
#export PATH="/edgehpc/dept/human_genetics/users/jjohn1/Discovery_Genewise_Merged/gwas2vcf/env/bin:$PATH"


library(dplyr)
library(data.table)
library(tools)
library(stringr)
library(glue)

options(scipen=999) 

eaf=fread("assocvariants.annotated_effectallelefrq.txt",sep="\t")
eaf$effectAlleleFreq<-as.numeric(eaf$effectAlleleFreq)

file="/edgehpc/dept/human_genetics/pQTL_deCODE_Ferkingstad_2021/gwas/5017_19_PRDX5_Peroxiredoxin_5.txt.gz"

n=0


files <- dir("/edgehpc/dept/human_genetics/pQTL_deCODE_Ferkingstad_2021/gwas", pattern = "\\.txt\\.gz$|\\.csv$", full.names = TRUE)


for ( file in files){
    print(glue(" now running {n}  {file}  "))

    df=fread(file,sep="\t")
    df=df[,c("Chrom","Pos","effectAllele","otherAllele","Beta","SE","minus_log10_pval","N","Name")]
    df$minus_log10_pval <- as.numeric(df$minus_log10_pval)
    df<- subset(df, minus_log10_pval >= 6)
    df$minus_log10_pval <- ifelse(df$minus_log10_pval < 200, df$minus_log10_pval, 199)
    df <- merge(df, eaf, by = "Name")
    df=df[,c("Chrom","Pos","Name","effectAllele","otherAllele","Beta","SE","N","minus_log10_pval","effectAlleleFreq")]
    df$minus_log10_pval <- 10^(-df$minus_log10_pval)
    df$Chrom <- gsub("chr", "", df$Chrom)

    outname<-basename(file)
    outname<-substr(outname, 1, nchar(outname) - 7)
    ID<-str_split(outname,"_")[[1]][3]

    fwrite(df, file = glue("{outname}.csv"), sep = "\t",row.names = FALSE)
    path=getwd()
    n=n+1

    r_command <- glue("python /edgehpc/dept/human_genetics/users/jjohn1/Discovery_Genewise_Merged/gwas2vcf/main.py \\
    --out {path}/{outname}.vcf \\
    --data {path}/{outname}.csv \\
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
    --dbsnp /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/GCF_000001405_cleaed.40_Msplited.vcf.gz \\
    --json params.txt --id {ID}")

    # Execute the R command using 'system'
    system(r_command)
    system('tabix -f -p vcf {path}/{outname}.vcf.gz')
}



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
    "build": "GRCh38"
    }
