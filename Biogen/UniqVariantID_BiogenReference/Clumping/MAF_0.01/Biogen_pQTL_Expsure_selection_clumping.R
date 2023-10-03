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
#/edgehpc/dept/human_genetics/users/jjohn1/Cis_Trans/Exposure/UKB_50KClimped/Reanalysis//mnt/depts/dept04/human_genetics/users/jjohn1/Cis_Trans/Exposure/UKB_50KClumped/Biogen/Reanalysis/Part0/All_Biogen_Exposure_TSS_chr_index_3k_Part0.csv"

files_df <- fread("All_Biogen_Exposure_TSS_chr_index_3k_Part0.csv")
tss<-files_df[,c("filenames","Symbol","gene", "Chromosome_edit","Biomarker_TSS_Start","Biomarker_TSS_End")]
tss$exposure <- gsub("^discovery_RS2_|\\.tsv\\.gz$", "", tss$filenames)

colnames(tss)<-c("Filename","Symbol","Gene name","Chromosome_edit","Biomarker_TSS_Start","Biomarker_TSS_End","exposure")

maf_cutoff<-0.01


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

no_cisvariant<-vector("character", length = 0)
no_transvariant<-vector("character", length = 0)

all_selected_cis_exposure_df<-data.frame()   #After pruning
all_selected_trans_exposure_df<-data.frame()  #After pruning
all_selected_trans_exposurenomhc_df<-data.frame()
all_significantPQTL_variant<-data.frame()  #Before pruning and qc

n<-1


for ( file in unique(tss$Filename) ){
    #file<-tss$Filename[1]
    print(glue("Now running {1} file :  {file}"))
    n<-n+1
    trans<-data.frame()
    cis<-data.frame()
    df<-data.frame()

    vcf_file<-glue("/edgehpc/dept/human_genetics/users/jjohn1/Pqtl_Summarystatistics/Biogen/VcfFiles/{substr(file, 1, nchar(file) - 7)}.vcf.gz")
    tryCatch({ vcf <- readVcf(vcf_file) }, error = function(e) {})
    tryCatch({df<-vcf_to_granges(vcf) %>% dplyr::as_tibble()}, error = function(e) {})
    tss2<-tss[tss$Filename==file][,c("Symbol","Chromosome_edit", "Biomarker_TSS_Start", "Biomarker_TSS_End","exposure")]

    if (nrow(df)>0) {
        df$LP <- as.numeric(df$LP)
        df["P"]<-10**-df$LP
        df<-df[df$P<=0.00000005,]
        df<-as.data.frame(df)

        collated_protein_dis <- merge(df, tss2, by.x = c("id"),by.y=c("Symbol"))
        # Calculate cis/trans status
        collated_protein_dis$Biomarker_TSS_Start <- as.numeric(collated_protein_dis$Biomarker_TSS_Start)
        collated_protein_dis$Biomarker_TSS_End <- as.numeric(collated_protein_dis$Biomarker_TSS_End)
        collated_protein_dis$start <- as.numeric(collated_protein_dis$start)

        collated_protein_dis$diff <- collated_protein_dis$Biomarker_TSS_Start - collated_protein_dis$start
        collated_protein_dis$diff2 <- collated_protein_dis$Biomarker_TSS_Start + collated_protein_dis$start

        collated_protein_dis$diff <- abs(collated_protein_dis$diff)
        collated_protein_dis$diff2 <- abs(collated_protein_dis$diff2)

        collated_protein_dis <- collated_protein_dis %>% mutate(Effect = case_when(
            seqnames == Chromosome_edit & (collated_protein_dis$diff < 1000000 | collated_protein_dis$diff2 < 1000000) ~ "CIS",
            seqnames == Chromosome_edit & (!collated_protein_dis$diff < 1000000 | collated_protein_dis$diff < 1000000) ~ "TRANS",
            !seqnames == Chromosome_edit ~ "TRANS"))

        all_significantPQTL_variant<-rbind.fill(all_significantPQTL_variant,collated_protein_dis)

        collated_protein_dis <- collated_protein_dis[(nchar(collated_protein_dis$REF) == 1 & nchar(collated_protein_dis$ALT) == 1),]
        collated_protein_dis$MAF <- ifelse(collated_protein_dis$AF > 0.5, 1-collated_protein_dis$AF, collated_protein_dis$AF)
        collated_protein_dis<-collated_protein_dis[!(((collated_protein_dis$REF %in% c("C", "G")) & (collated_protein_dis$ALT %in% c("C", "G"))) | ((collated_protein_dis$REF %in% c("A", "T")) & (collated_protein_dis$ALT %in% c("A", "T"))) & (collated_protein_dis$MAF > 0.42)), ]

        collated_protein_dis <- collated_protein_dis[collated_protein_dis$MAF >=maf_cutoff,]
        collated_protein_dis$ID <- paste(collated_protein_dis$seqnames, collated_protein_dis$start, collated_protein_dis$REF, collated_protein_dis$ALT, sep = "_")

        trans<-collated_protein_dis[collated_protein_dis$Effect=="TRANS",]
        cis<-collated_protein_dis[collated_protein_dis$Effect!="TRANS",]

        trans<-trans[,c("seqnames","start","ID","REF","ALT","ES","SE","P","AF","SI","SS","id","exposure")]
        cis<-cis[,c("seqnames","start","ID","REF","ALT","ES","SE","P","AF","SI","SS","id","exposure")] }

        if (nrow(cis)>0){
            selected_exposure_df<-data.frame()
            selected_exposure_df<-ld_pruning(cis)
            all_selected_cis_exposure_df<-rbind.fill(all_selected_cis_exposure_df,selected_exposure_df)
        }else {
        no_cisvariant<-c(no_cisvariant,file)
        cat(glue("\n\n{file} contains {nrow(cis)}\n\n"))}

        if (nrow(trans)>0){
            selected_exposure_df<-data.frame()
            selected_exposure_df<-ld_pruning(trans)
            all_selected_trans_exposure_df<-rbind.fill(all_selected_trans_exposure_df,selected_exposure_df)

            trans_nomhc <- trans[!(as.character(trans$seqnames) == "6" & trans$start >= 25000000 & trans$start <= 35000000), ]
            selected_exposure_omhc_df<-ld_pruning(trans_nomhc)
            all_selected_trans_exposurenomhc_df<-rbind.fill(all_selected_trans_exposurenomhc_df,selected_exposure_omhc_df)
        }else {
        no_transvariant<-c(no_transvariant,file)
            cat(glue("\n\n{file} contains {nrow(trans)}\n\n"))  }
    }

write.csv(all_selected_cis_exposure_df,glue("All_significannt_cis_exposure_AfterQC_LDclumping.csv"))
write.csv(all_selected_trans_exposure_df,glue("All_significannt_trans_exposure_AfterQC_LDclumping.csv"))
write.csv(all_selected_trans_exposurenomhc_df,glue("All_significannt_trans_exposure_AfterQC_LDclumping_MHCRemoval.csv"))
write.csv(all_significantPQTL_variant,glue("All_significannt_CisTransexposure_BeforeQC_LDclumping.csv"))

no_cisvariant_df <- data.frame(WithoutCisVariant = no_cisvariant)
write.csv(no_cisvariant_df,glue("GenesWithoutCisvariant.csv"))
no_transvariant_df <- data.frame(WithoutTransVariant = no_transvariant)
write.csv(no_transvariant_df,glue("GenesWithoutTransvariant.csv"))

