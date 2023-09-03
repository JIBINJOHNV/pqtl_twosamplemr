#.libPaths("/home/jjohn1/modulefiles/R4.3.0-foss-2021b")
#library(argparse)

.libPaths("/edgehpc/dept/human_genetics/jrobins_legacy/r_libs_4.1.0")
#.libPaths("/home/jjohn1/modulefiles/R4.3.0-foss-2021b")

suppressWarnings(suppressPackageStartupMessages({
    library(glue)
    library(gwasvcf)
    library(VariantAnnotation)
    library(dplyr)
    library(magrittr)
    library(VariantAnnotation)
    library(TwoSampleMR)
    library(MRInstruments)
    library(gwasglue)
    library(data.table)
    library(ieugwasr)
    library(plyr)
    library(mrpipeline)
    library(data.table)
}))

#parser <- ArgumentParser()
#parser$add_argument("--filename", help = "Provide read count file, it should only contain feature name and counts")
#args <- parser$parse_args()
#file_path <- glue("{args$filename}")

files<-c("ukb-d-30810_irnt_GRCh38_UniqID.vcf.gz")

# LD pruning
ld_pruning <- function(df2) {
        df3=data.frame()
        # Perform LD clumping and handle potential errors

    tryCatch({df3 <- suppressMessages(suppressWarnings({
                            ld_clump(dplyr::tibble(rsid = df2$ID, pval = df2$P, id = df2$id),
                            plink_bin = "/home/jjohn1/modulefiles/Softwares/Plink/plink",
                            bfile = "/edgehpc/dept/human_genetics/users/jjohn1/Referencefile_h19/ukb_imp_chr1_22_mac50_info07_b0_7_patched_bfiles_Hg38/ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles")  }))},
                        error = function(e) {# Handle the error
                        })

        if (nrow(df3)==0) {
            selected_df <- df2
        } else {
            df3 <- as.data.frame(df3)
            df3 <- df3[, c("rsid", "id")]  # Use "rsid" and "id" as column names in df3
            colnames(df3)[colnames(df3) == "id"] <- "id"
            colnames(df3)[colnames(df3) == "rsid"] <- "ID"
            selected_df <- merge(df2, df3, by = c("ID", "id"), all = FALSE)
        }
        return(selected_df) }

no_variant<-vector("character", length = 0)

all_selected_exposure_df<-data.frame()   #After pruning
all_selected_exposure_Nomhc_df<-data.frame()  #After pruning
all_significantPQTL_variant<-data.frame()  #Before pruning and qc

n<-1

file<-"ukb-d-30810_irnt_GRCh38_UniqID.vcf.gz"

for ( file in files) {
    #file<-tss$Filename[1]
    print(glue("Now running {1} file :  {file}"))
    n<-n+1
    df<-data.frame()
    vcf_file<-glue("/edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/Phosphate_ukb-d-30810_irnt/{file}")
    tryCatch({ vcf <- readVcf(vcf_file) }, error = function(e) {})
    tryCatch({df<-vcf_to_granges(vcf) %>% dplyr::as_tibble()}, error = function(e) {})

    if (nrow(df)>0) {
        df$LP <- as.numeric(df$LP)
        df["P"]<-10**-df$LP
        df<-df[df$P<=0.00000005,]
        collated_protein_dis<-as.data.frame(df)


        collated_protein_dis <- collated_protein_dis[(nchar(collated_protein_dis$REF) == 1 & nchar(collated_protein_dis$ALT) == 1),]
        collated_protein_dis$MAF <- ifelse(collated_protein_dis$AF > 0.5, 1-collated_protein_dis$AF, collated_protein_dis$AF)
        collated_protein_dis<-collated_protein_dis[!(((collated_protein_dis$REF %in% c("C", "G")) & (collated_protein_dis$ALT %in% c("C", "G"))) | ((collated_protein_dis$REF %in% c("A", "T")) & (collated_protein_dis$ALT %in% c("A", "T"))) & (collated_protein_dis$MAF > 0.42)), ]

        collated_protein_dis <- collated_protein_dis[collated_protein_dis$MAF >=0.001,]
        collated_protein_dis$ID <- paste(collated_protein_dis$seqnames, collated_protein_dis$start, collated_protein_dis$REF, collated_protein_dis$ALT, sep = "_")

        collated_protein_dis<-collated_protein_dis[,c("seqnames","start","ID","REF","ALT","ES","SE","P","AF","SI","SS","id")]
        collated_protein_dis$exposure<-collated_protein_dis$id

        all_significantPQTL_variant<-rbind.fill(all_significantPQTL_variant,collated_protein_dis)
        
        if (nrow(collated_protein_dis)>0){
            selected_exposure_df<-data.frame()
            selected_exposure_df<-ld_pruning(collated_protein_dis)
            all_selected_exposure_df<-rbind.fill(all_selected_exposure_df,selected_exposure_df)

            selected_exposure_Nomhc_df<-data.frame()
            selected_exposure_Nomhc_df <- collated_protein_dis[!(as.character(collated_protein_dis$seqnames) == "6" & collated_protein_dis$start >= 25000000 & collated_protein_dis$start <= 35000000), ]
            selected_exposure_Nomhc_df<-ld_pruning(selected_exposure_Nomhc_df)
            all_selected_exposure_Nomhc_df<-rbind.fill(all_selected_exposure_Nomhc_df,selected_exposure_Nomhc_df)

        }else {
        no_variant<-c(no_variant,file)
        cat(glue("\n\n{file} contains {nrow(cis)}\n\n"))}
    }
}


write.csv(all_selected_exposure_df,glue("All_significannt_exposure_AfterQC_LDclumping.csv"))
write.csv(all_selected_exposure_Nomhc_df,glue("All_significannt_exposure_AfterQC_LDclumping_MHCRemoval.csv"))
write.csv(all_significantPQTL_variant,glue("All_significannt_exposure_BeforeQC_LDclumping.csv"))

no_cisvariant_df <- data.frame(WithoutCisVariant = no_variant)
write.csv(no_cisvariant_df,glue("ExposureWithout_NOsignificantVariant.csv"))




#!/bin/bash

#SBATCH --ntasks=4
#SBATCH --time=00:01:00
#SBATCH --mem=120GB
#SBATCH --partition=cpu
#SBATCH --time=100:00:00

#source /home/jjohn1/modulefiles/anaconda3/bin/activate
#module load SAMtools/1.15-GCC-11.2.0
#module load PLINK/1.9b_6.21-x86_64
#module load R/4.1.3-foss-2021b
