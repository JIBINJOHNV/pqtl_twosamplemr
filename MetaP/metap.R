library(metap)
library(plyr)
library(glue)
library(ggplot2)
library(dplyr)



#matching_files <- list.files(path ="/Users/jibinjohn/Downloads/BiogenUKB_Decode_Pqtl_MR/Formeta/", pattern = "Decode_PGC3_SCZ", full.names = TRUE)


gwasnames <- c("ASD_PGC","BIP_PGC3_noukb","Depression_iPSYCH_2023","PGC3_SCZ","PGC_ADHD2022_iPSYCH_deCODE","PGC_AN2"  )

prefixes <- c( "CisExposure_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta",
               "TransExposureNoMHCUnique_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta",
                "TransExposureNoMHC_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta",
                "TransExposure_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta")


pqtltype1 <- "Decode"
pqtltype2 <- "Biogen"
exclude_columns <- c("Gene_Symbol", "outcome")

for (gwasname in gwasnames) {
  for (prefix in prefixes) {

    ## Decode PQTL
    decode <- read.csv(glue("{pqtltype1}_{gwasname}_{prefix}.csv"))
    decode <- setNames(decode, gsub("^exposure$", "Decodee.PQTL_exposure", colnames(decode)))

    ## Biogen PQTL
    biogen <- read.csv(glue("{pqtltype2}_{gwasname}_{prefix}.csv"))
    biogen <- setNames(biogen, gsub("^exposure$", "Biogen.PQTL_exposure", colnames(biogen)))

    ## Merge biogen and decode data frames
    decode_biogen <- merge(decode, biogen, by = c("Gene_Symbol", "outcome"), all = TRUE)

    #MRIVWtest_Pvalue
    MRIVWtest_Pvalue_cols<-colnames(decode_biogen[, grep("IVWDelta_Pvalue$", colnames(decode_biogen))])
    MRIVWtest_miss<-decode_biogen[!complete.cases(decode_biogen[MRIVWtest_Pvalue_cols]), ]
    MRIVWtest_nomiss<-decode_biogen[complete.cases(decode_biogen[MRIVWtest_Pvalue_cols]), ]

    MRIVWtest_df <- MRIVWtest_nomiss[, grep("IVWDelta_Pvalue$", colnames(decode_biogen))]
    MRIVWtest_df <- MRIVWtest_df %>%mutate_at(MRIVWtest_Pvalue_cols, as.numeric)

    metap_MRIVWtest<-vector()
    for(indx in 1:nrow(MRIVWtest_df)) {
            teest<-unlist(as.vector(MRIVWtest_df[indx,]))
            metap_MRIVWtest<-c(metap_MRIVWtest,allmetap(teest,method="sumlog")$p[[1]] )
            }
   
    MRIVWtest_nomiss[glue("{gwasname}_metap_MRIVWtest")]<-metap_MRIVWtest
    MRIVWtest_miss[glue("{gwasname}_metap_MRIVWtest")]<-"NA"
    decode_biogen <- rbind(MRIVWtest_nomiss,MRIVWtest_miss)

    out_name<-substr(prefix, 1, nchar(prefix) - 55)
    MRIVWtest_Beta_df <- MRIVWtest_nomiss[, grep("IVWDelta_Beta$", colnames(decode_biogen))]
    write.csv(MRIVWtest_df, glue("Biogen_Decode_pQTL_{gwasname}_{out_name}_Pvalue_MetapAnalysis.csv"), row.names = FALSE)
    system(glue("python3 altman_correlation.py -inputfilename Biogen_Decode_pQTL_{gwasname}_{out_name}_Pvalue_MetapAnalysis.csv -valutype pvalues ") )
    system(glue("rm Biogen_Decode_pQTL_{gwasname}_{out_name}_Pvalue_MetapAnalysis.csv") )

    write.csv(MRIVWtest_Beta_df, glue("Biogen_Decode_pQTL_{gwasname}_{out_name}_BetaValue_MetapAnalysis.csv"), row.names = FALSE)
    system(glue("python3 altman_correlation.py -inputfilename Biogen_Decode_pQTL_{gwasname}_{out_name}_BetaValue_MetapAnalysis.csv -valutype betavalues ") )
    system(glue("rm Biogen_Decode_pQTL_{gwasname}_{out_name}_BetaValue_MetapAnalysis.csv") )

    first_cols <- colnames(decode_biogen[, grep("IVWDelta_Pvalue$|Biogen_pval$", colnames(decode_biogen))])
    second_cols <- colnames(decode_biogen[, grep("_outcome$|_exposure$", colnames(decode_biogen))])
    third_cols <- colnames(decode_biogen[, grep("IVWDelta_Beta$", colnames(decode_biogen))])
    fourth_cols <- colnames(decode_biogen[, grep("IVWDelta_nSNPs$", colnames(decode_biogen))])
    fifth_cols <- colnames(decode_biogen[, grep("IVWDelta_SNPs$", colnames(decode_biogen))])
    sixth_cols <- colnames(decode_biogen[, grep("Heter.Stat$", colnames(decode_biogen))])
    seventh_cols <- colnames(decode_biogen[, grep("IVWDelta_SE$", colnames(decode_biogen))])

    new_order <- c(c("Gene_Symbol","outcome"), first_cols, c(glue("{gwasname}_metap_MRIVWtest")), second_cols,third_cols,fourth_cols,fifth_cols,sixth_cols,seventh_cols,
                    setdiff(colnames(decode_biogen), c("Gene_Symbol","outcome", first_cols,glue("{gwasname}_metap_MRIVWtest"), second_cols,third_cols,fourth_cols,fifth_cols,sixth_cols,seventh_cols)))

    decode_biogen_reordered_df <- decode_biogen[, new_order]

    if (grepl("CisExposure_" ,prefix )) {decode_biogen_reordered_df <- decode_biogen_reordered_df %>%rename_with(~ paste("CIS", ., sep = "_"), -c("Gene_Symbol", "outcome")) }
    if (grepl("TransExposure_British_",prefix )) {decode_biogen_reordered_df <- decode_biogen_reordered_df %>%rename_with(~ paste("TransExposure", ., sep = "_"), -c("Gene_Symbol", "outcome")) }
    if (grepl("TransExposureNoMHC_" ,prefix )) {decode_biogen_reordered_df <- decode_biogen_reordered_df %>%rename_with(~ paste("TransExposureNoMHC", ., sep = "_"), -c("Gene_Symbol", "outcome")) }
    if (grepl("TransExposureNoMHCUnique_" ,prefix )) {decode_biogen_reordered_df <- decode_biogen_reordered_df %>%rename_with(~ paste("TransExposureNoMHCUnique", ., sep = "_"), -c("Gene_Symbol", "outcome")) }

    write.csv(decode_biogen_reordered_df, glue("Biogen_Decode_pQTL_{gwasname}_{out_name}_MetapAnalysis.csv"), row.names = FALSE)


  }
}
