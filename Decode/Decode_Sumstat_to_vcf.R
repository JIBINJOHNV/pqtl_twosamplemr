library(dplyr)
library(data.table)
library(tools)
library(stringr)
library(glue)
suppressMessages(library(argparse))
options(scipen=999) 

parser <- ArgumentParser()


parser$add_argument("--filename", help = "Provide read count file, it should only contain feature name and counts")
args <- parser$parse_args()

file_path <- glue("{args$filename}")
files_df <- fread(file_path)




eaf=fread("/edgehpc/dept/human_genetics/users/jjohn1/Pqtl_Summarystatistics/DecodeIcelanders_35559_4907Protein/VcfFiles/assocvariants.annotated_effectallelefrq.txt",sep="\t")
eaf$effectAlleleFreq<-as.numeric(eaf$effectAlleleFreq)

n=0

for ( file in unique(files_df$Filename) ){
    print(glue("Now running {n} {file}"))
    df<-data.frame()
    outname<-substr(file, 1, nchar(file) - 7)
    ID<-str_split(outname,"_")[[1]][3]

    tryCatch({df=fread(glue("/edgehpc/dept/human_genetics/pQTL_deCODE_Ferkingstad_2021/gwas/redownload/{file}"),sep="\t")}, error = function(e) {})
    if (nrow(df)==0 ){tryCatch({df=fread(glue("/edgehpc/dept/human_genetics/pQTL_deCODE_Ferkingstad_2021/gwas/{file}"),sep="\t")}, error = function(e) {})}

    if ( nrow(df)>0  ){
            df=df[,c("Chrom","Pos","effectAllele","otherAllele","Beta","SE","minus_log10_pval","N","Name")]
            df$minus_log10_pval <- as.numeric(df$minus_log10_pval)
            df<- subset(df, minus_log10_pval >= 6)
            df$minus_log10_pval <- ifelse(df$minus_log10_pval < 200, df$minus_log10_pval, 199)
            df <- merge(df, eaf, by = "Name")
            df=df[,c("Chrom","Pos","Name","effectAllele","otherAllele","Beta","SE","N","minus_log10_pval","effectAlleleFreq")]
            df$minus_log10_pval <- 10^(-df$minus_log10_pval)
            df$Chrom <- gsub("chr", "", df$Chrom)

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
            system(glue('tabix -f -p vcf {path}/{outname}.vcf.gz'))

            }
}


#srun -n 12 --mem=30G --time=1:01:00  --pty bash
#module load Python/3.8.8-GCCcore-11.2.0
#module load libxml2/2.9.8-GCCcore-6.4.0
#module load OpenSSL/1.1
#source /edgehpc/dept/human_genetics/users/jjohn1/Discovery_Genewise_Merged/gwas2vcf/env/bin/activate
#export PATH="/edgehpc/dept/human_genetics/users/jjohn1/Discovery_Genewise_Merged/gwas2vcf/env/bin:$PATH"
#Rscript Decode_summary_to_vcf.R --filename /edgehpc/dept/human_genetics/users/jjohn1/Cis_Trans/1kg_clump/Decode/Decode_Proteins_TrnnscriptSTart_End_Hg38_Part1.csv
