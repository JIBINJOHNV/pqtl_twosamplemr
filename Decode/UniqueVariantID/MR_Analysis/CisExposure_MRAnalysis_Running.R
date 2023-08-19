.libPaths("/edgehpc/dept/human_genetics/jrobins_legacy/r_libs_4.1.0")

#.libPaths("/home/jjohn1/modulefiles/R4.1_modules")
#.libPaths("/home/jjohn1/modulefiles/anaconda3/envs/R43/lib/R/library")

#.libPaths("/home/jjohn1/modulefiles/R4.3.0-foss-2021b")
#.libPaths("/home/jjohn1/modulefiles/R4.1_modules")

suppressWarnings(suppressPackageStartupMessages({
    library(tidyr)
    library(glue)
    library(dplyr)
    library(plyr)
    library(data.table)
    library(magrittr)
    library(TwoSampleMR)
    library(ieugwasr)
    library(mrpipeline)
    library(gwasvcf)
    library(gwasglue)
    library(MRInstruments)
    library(VariantAnnotation)
}))

#install.packages("MendelianRandomization", lib="/home/jjohn1/modulefiles/mrpipeline/")

library(MendelianRandomization, lib.loc = '/home/jjohn1/modulefiles/mrpipeline')
#load("/home/jjohn1/modulefiles/R4.1_modules/rf.rdata")


exposure_file<-"All_significannt_cis_exposure_AfterQC_LDclumping.csv"
prefix="cis_exposure"

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
                eaf_col = "AF",
                effect_allele_col = "ALT",other_allele_col = "REF", # nolint: line_length_linter.
                pval_col = "P",
                samplesize_col = "SS",
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


###MR Analysis
perform_mr_analysis <- function(dat) {
        mr_res <- data.frame()
        mr_het <- data.frame()
        mr_pleo <- data.frame()
        mr_single_wald_ratio <- data.frame()
        mr_single_meta_fixed <- data.frame()
        mr_direction <- data.frame()

        tryCatch({
            mr_res <- mr(dat, method_list = methodlist_1)
            mr_res <- generate_odds_ratios(mr_res)
            mr_res_df <<- rbind.fill(mr_res_df, mr_res)
        }, error = function(e) {
            # Handle the error or perform any desired action
        })

        if (nrow(mr_res)==0) {
        tryCatch({
            mr_res <- mr(dat, method_list = methodlist_2)
            mr_res <- generate_odds_ratios(mr_res)
            mr_res_df <<- rbind.fill(mr_res_df, mr_res)
        }, error = function(e) {
            # Handle the error or perform any desired action
        })}

        if (nrow(mr_res)==0) {
        tryCatch({
            mr_res <- mr(dat, method_list = methodlist_3)
            mr_res <- generate_odds_ratios(mr_res)
            mr_res_df <<- rbind.fill(mr_res_df, mr_res)
        }, error = function(e) {
            # Handle the error or perform any desired action
        })}

        tryCatch({
            mr_het <- mr_heterogeneity(dat, method_list = heterogenoty_list_1)
        }, error = function(e) {})

        if (nrow(mr_het)==0) {
        tryCatch({
            mr_het <- mr_heterogeneity(dat, method_list = heterogenoty_list_2)
        }, error = function(e) {
            # Handle the error or perform any desired action
        })}

        if (nrow(mr_het)>0) {
            mr_het_df <<- rbind.fill(mr_het_df, mr_het)
        }

       tryCatch({
            mr_pleo <- mr_pleiotropy_test(dat)
        }, error = function(e) {})

        if (nrow((mr_pleo)>0)) {
            mr_pleo_df <<- rbind.fill(mr_pleo_df, mr_pleo)
        }

        tryCatch({
            mr_direction <- directionality_test(dat)
        }, error = function(e) {})

        if (nrow(mr_direction)>0) {
            direction_df <<- rbind.fill(direction_df, mr_direction)
        }}



##Single Variant annalysis
perform_single_variant_mr_analysis <- function(dat) {
        mr_single_wald_ratio <- NULL
        mr_single_meta_fixed <- NULL

        tryCatch({
            mr_single_wald_ratio <- mr_singlesnp(dat, parameters = default_parameters(), single_method = "mr_wald_ratio")
            }, error = function(e) {})

        if ( !is.null(mr_single_wald_ratio) || nrow(mr_res) != 0) {wald_df <<- rbind.fill(wald_df, mr_single_wald_ratio)}

        tryCatch({
            mr_single_meta_fixed <- mr_singlesnp(dat, parameters = default_parameters(), single_method = "mr_meta_fixed")
        }, error = function(e) {})

        if ( !is.null(mr_single_meta_fixed) || nrow(mr_res) != 0  ) { meta_fixed_df <<- rbind.fill(meta_fixed_df, mr_single_meta_fixed)}
            }


##Create empty data frame
no_outcome_dat_df<-data.frame() ##No varaints or proxy was present in the outcome data
harmonised_df<-data.frame()
mr_res_df <- data.frame()
mr_het_df <- data.frame()
mr_pleo_df <- data.frame()
direction_df <- data.frame()
meta_fixed_df <- data.frame()
wald_df <- data.frame()

mendelianpipeline_mrresult<-data.frame()
MendelianRandomization_df<-data.frame()
mr_presso_df<-data.frame()
mr_ivw_delta_df<-data.frame()

###Exposure file
exposure_df <- fread(exposure_file)
if ("V1" %in% colnames(exposure_df)) { exposure_df <- exposure_df[, -"V1", with = FALSE] }

if ("X" %in% unique(exposure_df$seqnames)) { exposure_df<-exposure_df[exposure_df$seqnames!="X"] }



for ( Exposure in sort(unique(exposure_df$exposure))  ) {

        formated_outcomedf <- data.frame()
        formated_selected_exposure_df <- data.frame()
        formated_outcomedf <- data.frame()
        dat <- data.frame()
        df2 <- data.frame()
        MRAllObject_all <- data.frame()
        results <- data.frame()
        mr_res <- data.frame()
        dat2 <-  data.frame()

        tryCatch({
            exposure_selected_df <- exposure_df[exposure_df$exposure == as.character(Exposure), ]
        }, error = function(e) {})

        if (nrow(exposure_selected_df) > 0) {
            formated_selected_exposure_df <- format_exposure(exposure_selected_df)
        }

        if (nrow(formated_selected_exposure_df) > 0) {
            formated_selected_exposure_df <- calc_f_stat(formated_selected_exposure_df, f_cutoff = 0)
            formated_selected_exposure_df$SNP <- toupper(formated_selected_exposure_df$SNP)

        }

        tryCatch({
            if (nrow(formated_selected_exposure_df) > 0) {
                formated_selected_exposure_df$exposure<-Exposure
                formated_outcomedf <-  format_outcome(outcomevcf, outcomevcf_index,formated_selected_exposure_df, dbfile)
                formated_outcomedf$SNP <- toupper(formated_outcomedf$SNP)
            }
        }, error = function(e) {})

        if (nrow(formated_outcomedf) > 0) {
            dat <- harmonise_data(exposure_dat = formated_selected_exposure_df, outcome_dat = formated_outcomedf,action = 1)
        } else {
            no_outcome_dat_df <- rbind.fill(no_outcome_dat_df, formated_selected_exposure_df)
        }

        #Two Sample MR test https://mrcieu.github.io/TwoSampleMR/
        if (nrow(dat) > 0) {
            perform_mr_analysis(dat)
            harmonised_df <<- rbind.fill(harmonised_df, dat)
            perform_single_variant_mr_analysis(dat)
        }

        #Biogen mrpipeline; https://github.com/jwr-git/mrpipeline
        if (nrow(dat) > 0) {
            mr_res <- do_mr(dat, f_cutoff = 0, all_wr = TRUE, verbose = TRUE)
            mendelianpipeline_mrresult <<- rbind.fill(mendelianpipeline_mrresult, mr_res)
        }



        #MendelianRandomization  https://github.com/cran/MendelianRandomization Version 0.8
        if (nrow(dat) > 0) {
            M_Randomization_df<-data.frame()
            dat2 <- dat_to_MRInput(dat)
            tryCatch({
            MRAllObject_all <- MendelianRandomization::mr_allmethods(dat2[[1]], method = "all")
            }, error = function(e) {})

            tryCatch({
            M_Randomization_df <- as.data.frame(MRAllObject_all@Values)
            }, error = function(e) {})


            #Some situation all meethods not work , however ivw will work
            if (nrow(M_Randomization_df) == 0) {
            tryCatch({
            MRAllObject_all <- MendelianRandomization::mr_allmethods(dat2[[1]], method = "ivw")
            }, error = function(e) {})

            tryCatch({
            M_Randomization_df <- as.data.frame(MRAllObject_all@Values)
            }, error = function(e) {})

            }

            if (nrow(M_Randomization_df) > 0) {
            M_Randomization_df$outcome <- dat$outcome[1]
            M_Randomization_df$exposure <- dat$exposure[1]
            M_Randomization_df$SNPs <- paste(unique(dat$SNP), collapse = ",")
            MendelianRandomization_df <<- rbind.fill(MendelianRandomization_df, M_Randomization_df)
            }
        }



         ##presso

        if (nrow(dat) > 0) {
            mr_presso <- data.frame()
            presso <- data.frame()
            tryCatch({
            mr_presso <- run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
            }, error = function(e) {})

            tryCatch({
            presso <- mr_presso[[1]]$`Main MR results`
            }, error = function(e) {})

            if (nrow(presso) > 0) {
            presso$RSSobs <- mr_presso[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs
            presso$GlobalTest_Pvalue <- mr_presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
            presso$exposure <- attributes(mr_presso)$exposure
            presso$outcome <- attributes(mr_presso)$outcome
            mr_presso_df <<- rbind.fill(mr_presso_df, presso)
            }
        }

        if (nrow(dat) > 0) {
            t <- data.frame()
            tryCatch({
            t <- MendelianRandomization::mr_ivw(dat2[[1]], weights = "delta", distribution = "normal")
            }, error = function(e) {})

            if (typeof(t) == "S4") {
            mr_ivw_delta_df <<- rbind.fill(mr_ivw_delta_df, data.frame(
                Model = t@Model,
                Outcome = t@Outcome,
                Penalized = t@Penalized,
                Estimate = t@Estimate,
                CILower = t@CILower,
                Alpha = t@Alpha,
                SNPs = t@SNPs,
                Heter.Stat = paste(unique(t@Heter.Stat),collapse = ","),
                Exposure = t@Exposure,
                Robust = t@Robust,
                Correlation = t@Correlation,
                StdError = t@StdError,
                CIUpper = t@CIUpper,
                Pvalue = t@Pvalue,
                RSE = t@RSE
            ))
            }
        }
        }


write.csv(no_outcome_dat_df, glue('{prefix}_Novariants_In_Outcome_GWAS.csv'), row.names = FALSE)
write.csv(harmonised_df, glue('{prefix}_Harmonised_Exposure_Outcome.csv'),row.names = FALSE)
write.csv(mr_res_df, glue('{prefix}_TwoSampleMR_Analysis_Multiple_MR_Test.csv'),row.names = FALSE)
write.csv(mr_het_df, glue('{prefix}_TwoSampleMR_Analysis_HeterogenityTest.csv'),row.names = FALSE)
write.csv(mr_pleo_df, glue('{prefix}_TwoSampleMR_Analysis_Hpleiotropy_Test.csv'),row.names = FALSE)
write.csv(direction_df, glue('{prefix}_TwoSampleMR_Analysis_Directionality_Test.csv'),row.names = FALSE)
write.csv(meta_fixed_df, glue('{prefix}_TwoSampleMR_Analysis_SingleVariantMetafixed_Test.csv'),row.names = FALSE)
write.csv(wald_df, glue('{prefix}_TwoSampleMR_Analysis_SingleVariantWald_Test.csv'),row.names = FALSE)
write.csv(mendelianpipeline_mrresult, glue('{prefix}_MendelianPipelineTest.csv'),row.names = FALSE)
write.csv(MendelianRandomization_df, glue('{prefix}_MendelianRandomization_AllTest.csv'),row.names = FALSE)
write.csv(mr_presso_df, glue('{prefix}_MendelianRandomization_Presso_Test.csv'),row.names = FALSE)
write.csv(mr_ivw_delta_df, glue('{prefix}_MendelianRandomization_IVW_Delta_Test.csv'),row.names = FALSE)
