library(metap)
library(plyr)
library(glue)

gwasnames<-c('Depression_iPSYCH_2023','BIP_PGC3_noukb')
prefixes<-c("CisExposure_CompleteMR_AnalysisResults_WithselectedColumns_ForMeta.csv","TransExposureNoMHC_CompleteMR_AnalysisResults_WithselectedColumns_ForMeta.csv",
            "TransExposureNoMHCUnique_CompleteMR_AnalysisResults_WithselectedColumns_ForMeta.csv","TransExposure_CompleteMR_AnalysisResults_WithselectedColumns_ForMeta.csv")


pqtltype1<-"Decode"
pqtltype2<-"Biogen"
gwasname<- 'Depression_iPSYCH_2023'     
prefix<-"CisExposure_CompleteMR_AnalysisResults_WithselectedColumns_ForMeta"

exclude_columns <- c("Gene_Symbol","outcome")

##Decode PQTL 
decode <- read.csv(glue("{pqtltype1}_{gwasname}_{prefix}.csv"))
col_names <- colnames(decode)
columns_to_rename <- col_names[!(colnames(decode) %in% exclude_columns)]
new_col_names <- ifelse(columns_to_rename != "", paste0("DecodePqtl:", columns_to_rename), columns_to_rename)
colnames(decode)[!(col_names %in% exclude_columns)] <- new_col_names

##Biogen PQTL
biogen=read.csv(glue("{pqtltype2}_{gwasname}_{prefix}.csv"))
col_names <- colnames(biogen)
columns_to_rename <- col_names[!(colnames(biogen) %in% exclude_columns)]
new_col_names <- ifelse(columns_to_rename != "", paste0("BiogenPqtl:", columns_to_rename), columns_to_rename)
colnames(biogen)[!(col_names %in% exclude_columns)] <- new_col_names


##Merge biogen and decode dataframe
decode_biogen <- merge(decode, biogen, by = c("Gene_Symbol","outcome" ),all = TRUE)



#MRIVWtest_Pvalue
MRIVWtest_Pvalue_cols<-colnames(decode_biogen[, grep("IVWDelta_Pvalue$", colnames(decode_biogen))])
MRIVWtest_miss<-decode_biogen[!complete.cases(decode_biogen[MRIVWtest_Pvalue_cols]), ]
MRIVWtest_nomiss<-decode_biogen[complete.cases(decode_biogen[MRIVWtest_Pvalue_cols]), ]

MRIVWtest_df <- MRIVWtest_nomiss[, grep("IVWDelta_Pvalue$", colnames(decode_biogen))]

metap_MRIVWtest<-vector()

for(indx in 1:nrow(MRIVWtest_df)) {
  teest<-unlist(as.vector(MRIVWtest_df[indx,]))
  metap_MRIVWtest<-c(metap_MRIVWtest,allmetap(teest,method="sumlog")$p[[1]] )
}

MRIVWtest_nomiss$metap_MRIVWtest<-metap_MRIVWtest
MRIVWtest_miss$metap_MRIVWtest<-"NA"
decode_biogen <- rbind(MRIVWtest_nomiss,MRIVWtest_miss)


#MR_Pipeline_pval
MR_Pipeline_pval_cols<-colnames(decode_biogen[, grep("Biogen_pval$", colnames(decode_biogen))])
MR_Pipeline_miss<-decode_biogen[!complete.cases(decode_biogen[MR_Pipeline_pval_cols]), ]
MR_Pipeline_nomiss<-decode_biogen[complete.cases(decode_biogen[MR_Pipeline_pval_cols]), ]

MR_Pipeline_df <- MR_Pipeline_nomiss[, grep("Biogen_pval$", colnames(decode_biogen))]
metap_MR_Pipeline<-vector()

for(indx in 1:nrow(MR_Pipeline_df)) {
  teest<-unlist(as.vector(MR_Pipeline_df[indx,]))
  metap_MR_Pipeline<-c(metap_MR_Pipeline,allmetap(teest,method="sumlog")$p[[1]] )
}


MR_Pipeline_nomiss$metap_BiogenMR_Pipeline<-metap_MR_Pipeline
MR_Pipeline_miss$metap_BiogenMR_Pipeline<-"NA"
decode_biogen <- rbind(MR_Pipeline_nomiss,MR_Pipeline_miss)


##Reordering the columns
first_cols <- colnames(decode_biogen[, grep("IVWDelta_Pvalue$|Biogen_pval$", colnames(decode_biogen))])
second_cols <- colnames(decode_biogen[, grep("_Outcome$|_Exposure$", colnames(decode_biogen))])

new_order <- c(
  c("Gene_Symbol"), first_cols, c("metap_BiogenMR_Pipeline", "metap_MRIVWtest"), second_cols,
  setdiff(colnames(decode_biogen), c("Gene_Symbol", first_cols, "metap_BiogenMR_Pipeline", "metap_MRIVWtest", second_cols))
)


decode_biogen_reordered_df <- decode_biogen[, new_order]
write.csv(decode_biogen_reordered_df, glue("Biogen_Decode_pQTL_{gwasname}_{prefix}.csv"), row.names = FALSE)




library(metap)
library(plyr)
library(glue)

gwasnames <- c('Depression_iPSYCH_2023', 'BIP_PGC3_noukb')
prefixes <- c("CisExposure_CompleteMR_AnalysisResults_WithselectedColumns_ForMeta", "TransExposureNoMHC_CompleteMR_AnalysisResults_WithselectedColumns_ForMeta",
              "TransExposureNoMHCUnique_CompleteMR_AnalysisResults_WithselectedColumns_ForMeta", "TransExposure_CompleteMR_AnalysisResults_WithselectedColumns_ForMeta")

pqtltype1 <- "Decode"
pqtltype2 <- "Biogen"
exclude_columns <- c("Gene_Symbol", "outcome")

for (gwasname in gwasnames) {
  for (prefix in prefixes) {
    ## Decode PQTL
    decode <- read.csv(glue("{pqtltype1}_{gwasname}_{prefix}.csv"))
    col_names <- colnames(decode)
    columns_to_rename <- col_names[!(col_names %in% exclude_columns)]
    new_col_names <- ifelse(columns_to_rename != "", paste0("DecodePqtl:", columns_to_rename), columns_to_rename)
    colnames(decode)[!(col_names %in% exclude_columns)] <- new_col_names

    ## Biogen PQTL
    biogen <- read.csv(glue("{pqtltype2}_{gwasname}_{prefix}.csv"))
    col_names <- colnames(biogen)
    columns_to_rename <- col_names[!(colnames(biogen) %in% exclude_columns)]
    new_col_names <- ifelse(columns_to_rename != "", paste0("BiogenPqtl:", columns_to_rename), columns_to_rename)
    colnames(biogen)[!(col_names %in% exclude_columns)] <- new_col_names

    ## Merge biogen and decode data frames
    decode_biogen <- merge(decode, biogen, by = c("Gene_Symbol", "outcome"), all = TRUE)

    # Rest of your code for meta-analysis, reordering columns, and writing to CSV

    #MRIVWtest_Pvalue
    MRIVWtest_Pvalue_cols<-colnames(decode_biogen[, grep("IVWDelta_Pvalue$", colnames(decode_biogen))])
    MRIVWtest_miss<-decode_biogen[!complete.cases(decode_biogen[MRIVWtest_Pvalue_cols]), ]
    MRIVWtest_nomiss<-decode_biogen[complete.cases(decode_biogen[MRIVWtest_Pvalue_cols]), ]

    MRIVWtest_df <- MRIVWtest_nomiss[, grep("IVWDelta_Pvalue$", colnames(decode_biogen))]

    metap_MRIVWtest<-vector()

    for(indx in 1:nrow(MRIVWtest_df)) {
    teest<-unlist(as.vector(MRIVWtest_df[indx,]))
    metap_MRIVWtest<-c(metap_MRIVWtest,allmetap(teest,method="sumlog")$p[[1]] )
    }

    MRIVWtest_nomiss$metap_MRIVWtest<-metap_MRIVWtest
    MRIVWtest_miss$metap_MRIVWtest<-"NA"
    decode_biogen <- rbind(MRIVWtest_nomiss,MRIVWtest_miss)


    #MR_Pipeline_pval
    MR_Pipeline_pval_cols<-colnames(decode_biogen[, grep("Biogen_pval$", colnames(decode_biogen))])
    MR_Pipeline_miss<-decode_biogen[!complete.cases(decode_biogen[MR_Pipeline_pval_cols]), ]
    MR_Pipeline_nomiss<-decode_biogen[complete.cases(decode_biogen[MR_Pipeline_pval_cols]), ]

    MR_Pipeline_df <- MR_Pipeline_nomiss[, grep("Biogen_pval$", colnames(decode_biogen))]
    metap_MR_Pipeline<-vector()

    for(indx in 1:nrow(MR_Pipeline_df)) {
    teest<-unlist(as.vector(MR_Pipeline_df[indx,]))
    metap_MR_Pipeline<-c(metap_MR_Pipeline,allmetap(teest,method="sumlog")$p[[1]] )
    }


    MR_Pipeline_nomiss$metap_BiogenMR_Pipeline<-metap_MR_Pipeline
    MR_Pipeline_miss$metap_BiogenMR_Pipeline<-"NA"
    decode_biogen <- rbind(MR_Pipeline_nomiss,MR_Pipeline_miss)
    ##Reordering the columns
    first_cols <- colnames(decode_biogen[, grep("IVWDelta_Pvalue$|Biogen_pval$", colnames(decode_biogen))])
    second_cols <- colnames(decode_biogen[, grep("_Outcome$|_Exposure$", colnames(decode_biogen))])

    new_order <- c(
    c("Gene_Symbol"), first_cols, c("metap_BiogenMR_Pipeline", "metap_MRIVWtest"), second_cols,
    setdiff(colnames(decode_biogen), c("Gene_Symbol", first_cols, "metap_BiogenMR_Pipeline", "metap_MRIVWtest", second_cols))
    )

    decode_biogen_reordered_df <- decode_biogen[, new_order]
    write.csv(decode_biogen_reordered_df, glue("Biogen_Decode_pQTL_{gwasname}_{prefix}.csv"), row.names = FALSE)

    }}
