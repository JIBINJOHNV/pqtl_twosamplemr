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
library(MendelianRandomization, lib.loc = '/home/jjohn1/modulefiles/mrpipeline')

#parser <- ArgumentParser()
#parser$add_argument("--filename", help = "Provide read count file, it should only contain feature name and counts")
#args <- parser$parse_args()
#file_path <- glue("{args$filename}")
#/edgehpc/dept/human_genetics/users/jjohn1/Cis_Trans/Exposure/UKB_50KClimped/Reanalysis//mnt/depts/dept04/human_genetics/users/jjohn1/Cis_Trans/Exposure/UKB_50KClumped/Biogen/Reanalysis/Part0/All_Biogen_Exposure_TSS_chr_index_3k_Part0.csv"

files_df <- fread("All_significannt_CisTransexposure_BeforeQC_LDclumping.csv")
files_df$exposure <- as.character(files_df$exposure)





prefix="Biogen_PGC3_SCZ_CisExposure_test"

set_bcftools('/home/jjohn1/modulefiles/anaconda3/bin/bcftools')
set_plink('/home/jjohn1/modulefiles/Softwares/Plink/plink')

outcomevcf="/edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_Processed_Outcome_GWAS_UniqIds/PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf.gz"
outcomevcf_index="/edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_Processed_Outcome_GWAS_UniqIds/PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf.index"
dbfile="/edgehpc/dept/human_genetics/users/jjohn1/Referencefile_h19/ukb_imp_chr1_22_mac50_info07_b0_7_patched_bfiles_Hg38/ukb_imp_chr1_22XXY_mac50_info07_b0_7_patched_bfiles"

methodlist_1<-c("mr_wald_ratio","mr_two_sample_ml","mr_egger_regression","mr_egger_regression_bootstrap","mr_simple_median","mr_weighted_median","mr_penalised_weighted_median","mr_ivw","mr_ivw_radial","mr_ivw_mre","mr_ivw_fe","mr_simple_mode","mr_weighted_mode","mr_weighted_mode_nome","mr_simple_mode_nome","mr_sign","mr_uwr")
methodlist_2<-c("mr_wald_ratio","mr_two_sample_ml","mr_egger_regression","mr_egger_regression_bootstrap","mr_simple_median","mr_weighted_median","mr_penalised_weighted_median","mr_ivw","mr_ivw_mre","mr_ivw_fe","mr_simple_mode","mr_weighted_mode","mr_weighted_mode_nome","mr_simple_mode_nome","mr_sign","mr_uwr")
methodlist_3<-c("mr_wald_ratio","mr_two_sample_ml","mr_egger_regression","mr_egger_regression_bootstrap","mr_simple_median","mr_weighted_median","mr_penalised_weighted_median","mr_ivw")

heterogenoty_list_1<-c("mr_two_sample_ml","mr_egger_regression","mr_ivw","mr_ivw_radial","mr_uwr")
heterogenoty_list_2<-c("mr_two_sample_ml","mr_egger_regression","mr_ivw","mr_uwr")
#create_rsidx_index_from_vcf(outcomevcf, outcomevcf_index)



##Format Exposure
format_exposure <- function(selecteed_df) {
            exp_df <- format_data(
                selecteed_df,
                type = "exposure",snps = NULL,header = TRUE,
                phenotype_col = "exposure",snp_col = "ID",
                beta_col = "ES",se_col = "SE",
                eaf_col = "AF",effect_allele_col = "ALT",other_allele_col = "REF", # nolint: line_length_linter.
                pval_col = "P",samplesize_col = "SS",
                gene_col = "id",min_pval = 1e-800,
                chr_col = "seqnames",pos_col = "start",
                log_pval = FALSE
            )
            return(exp_df) }


#format Outcome data
format_outcome<-function(outcome_file,outcome_fileindex,exposure_df,ldfile){
                        vcf <- query_gwas(outcome_file,
                                        rsid=exposure_df$SNP,
                                        rsidx=outcome_fileindex,
                                        proxies="yes",
                                        bfile=ldfile,
                                        tag_r2=0.8)
                        outcome_df<-gwasvcf_to_TwoSampleMR(vcf, type = "outcome")
                        return(outcome_df)}


#Merge the results
merge_ivwdelta<-function(t_delta,cis_mr_ivw_delta_df,exposure_){
            if (typeof(t_delta) == "S4") {
                    cis_mr_ivw_delta_df <<- rbind.fill(cis_mr_ivw_delta_df, data.frame(
                        Model = t_delta@Model,Outcome = t_delta@Outcome,Penalized = t_delta@Penalized,Estimate = t_delta@Estimate,CILower = t_delta@CILower,Alpha = t_delta@Alpha,
                        SNPs = t_delta@SNPs,Heter.Stat = paste(unique(t_delta@Heter.Stat),collapse = ","),Exposure =exposure_,Robust = t_delta@Robust,StdError = t_delta@StdError,
                        CIUpper = t_delta@CIUpper,Pvalue = t_delta@Pvalue,RSE = t_delta@RSE ))} }


cis_mr_ivw_delta_df<-data.frame()
trans_mr_ivw_delta_df<-data.frame()


for (exposure_ in unique(files_df$exposure)) {
    collated_protein_dis <- files_df[files_df$exposure == exposure_]

    collated_protein_dis <- collated_protein_dis[(nchar(collated_protein_dis$REF) == 1 & nchar(collated_protein_dis$ALT) == 1), ]
    collated_protein_dis$MAF <- ifelse(collated_protein_dis$AF > 0.5, 1 - collated_protein_dis$AF, collated_protein_dis$AF)
    collated_protein_dis <- collated_protein_dis[!(((collated_protein_dis$REF %in% c("C", "G")) & (collated_protein_dis$ALT %in% c("C", "G"))) | ((collated_protein_dis$REF %in% c("A", "T")) & (collated_protein_dis$ALT %in% c("A", "T"))) & (collated_protein_dis$MAF > 0.42)), ]
    collated_protein_dis <- collated_protein_dis[collated_protein_dis$MAF >= 0.001, ]
    collated_protein_dis$ID <- paste(collated_protein_dis$seqnames, collated_protein_dis$start, collated_protein_dis$REF, collated_protein_dis$ALT, sep = "_")

    trans <- collated_protein_dis[collated_protein_dis$Effect == "TRANS", ]
    cis <- collated_protein_dis[collated_protein_dis$Effect != "TRANS", ]
    data_frames <- list(trans, cis)


    for (df in data_frames) {
        if (nrow(df) > 0) {
            formated_selected_exposure_df <- format_exposure(df)
        } else {
            no_cisvariant <- c(no_cisvariant, file)
            cat(glue("\n\n{file} contains {nrow(df)}\n\n"))
        }

        if (nrow(formated_selected_exposure_df) > 0) {
            formated_selected_exposure_df <- calc_f_stat(formated_selected_exposure_df, f_cutoff = 0)
            formated_selected_exposure_df$SNP <- toupper(formated_selected_exposure_df$SNP)
        }

        tryCatch({
            if (nrow(formated_selected_exposure_df) > 0) {
                formated_selected_exposure_df$exposure <- exposure_
                formated_outcomedf <- format_outcome(outcomevcf, outcomevcf_index, formated_selected_exposure_df, dbfile)
                formated_outcomedf$SNP <- toupper(formated_outcomedf$SNP)
            }
        }, error = function(e) {})

        if (nrow(formated_outcomedf) > 0) {
            dat <- harmonise_data(exposure_dat = formated_selected_exposure_df, outcome_dat = formated_outcomedf, action = 1)
        } else {
            no_outcome_dat_df <- rbind.fill(no_outcome_dat_df, formated_selected_exposure_df)
        }

        if (nrow(dat) > 1) {
            ldmatrix <- data.frame()
            MRInputObject.cor <- data.frame()
            t_delta <- data.frame()
            tryCatch({
                ldmatrix <- ld_matrix(dat$SNP, with_alleles = TRUE, bfile = dbfile, plink_bin = '/home/jjohn1/modulefiles/Softwares/Plink/plink')
                MRInputObject.cor <- mr_input(bx = dat$`beta.exposure`, bxse = dat$`se.exposure`, by = dat$`beta.outcome`, byse = dat$`se.outcome`, corr = ldmatrix)
                t_delta <- MendelianRandomization::mr_ivw(MRInputObject.cor, weights = "delta", distribution = "normal", correl = TRUE)

                if (deparse(substitute(cis)) == "cis") {
                    result <- merge_ivwdelta(t_delta, cis_mr_ivw_delta_df,exposure_)
                } else {
                    result <- merge_ivwdelta(t_delta, trans_mr_ivw_delta_df,exposure_)
                }
            }, error = function(e) {
                # Handle the error here if needed
            })
        }

        if (nrow(dat) == 1) {
            ldmatrix <- data.frame()
            MRInputObject.cor <- data.frame()
            t_delta <- data.frame()
            tryCatch({
                ldmatrix <- ld_matrix(dat$SNP, with_alleles = TRUE, bfile = dbfile, plink_bin = '/home/jjohn1/modulefiles/Softwares/Plink/plink')
                MRInputObject.cor <- mr_input(bx = dat$`beta.exposure`, bxse = dat$`se.exposure`, by = dat$`beta.outcome`, byse = dat$`se.outcome`, corr = ldmatrix)
                t_delta <- MendelianRandomization::mr_ivw(MRInputObject.cor, weights = "delta", distribution = "normal", correl = TRUE)

                if (deparse(substitute(cis)) == "cis") {
                    result <- merge_ivwdelta(t_delta, cis_mr_ivw_delta_df,exposure_)
                } else {
                    result <- merge_ivwdelta(t_delta, trans_mr_ivw_delta_df,exposure_)
                }
            }, error = function(e) {
                # Handle the error here if needed
            })
        }
    }
}
