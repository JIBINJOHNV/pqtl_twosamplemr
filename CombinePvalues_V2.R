
library(metap)
library(plyr)

prefix="trans_exposure_noMHC"
decode=read.csv("Decode_rans_exposure_noMHC_CompleteMR_AnalysisResults_ForMtaP.csv")
colnames(decode)[colnames(decode) == "Outcome"] <- "Decode_trans_exposure_noMHC_Unique_Outcome"
colnames(decode)[colnames(decode) == "Exposure"] <- "Decode_trans_exposure_noMHC_Unique_Exposure"


biogen=read.csv("Biogen_trans_exposure_NoMHC_CompleteMR_AnalysisResults_ForMtaP.csv")
colnames(biogen)<-paste("Biogen_trans_exposure_noMHC", colnames(biogen), sep = "_")
colnames(biogen)[colnames(biogen) == "Biogen_trans_exposure_noMHC_Gene_Symbol"] <- "Gene_Symbol"

#colnames(biogen)[colnames(biogen) == "Outcome"] <- "Biogen_trans_exposure_noMHC_Outcome"
#colnames(biogen)[colnames(biogen) == "Exposure"] <- "Biogen_trans_exposure_noMHC_Exposure"


decode_biogen <- merge(decode, biogen, by = "Gene_Symbol",all = TRUE)



#MRIVWtest_Pvalue
MRIVWtest_Pvalue_cols<-colnames(decode_biogen[, grep("MRIVWtest_Pvalue$", colnames(decode_biogen))])

MRIVWtest_miss<-decode_biogen[!complete.cases(decode_biogen[MRIVWtest_Pvalue_cols]), ]
MRIVWtest_nomiss<-decode_biogen[complete.cases(decode_biogen[MRIVWtest_Pvalue_cols]), ]


MRIVWtest_df <- MRIVWtest_nomiss[, grep("MRIVWtest_Pvalue$", colnames(decode_biogen))]

metap_MRIVWtest<-vector()

for(indx in 1:nrow(MRIVWtest_df)) {
  teest<-unlist(as.vector(MRIVWtest_df[indx,]))
  metap_MRIVWtest<-c(metap_MRIVWtest,allmetap(teest,method="sumlog")$p[[1]] )
}

MRIVWtest_nomiss$metap_MRIVWtest<-metap_MRIVWtest
MRIVWtest_miss$metap_MRIVWtest<-"NA"
decode_biogen <- rbind(MRIVWtest_nomiss,MRIVWtest_miss)



#MR_Pipeline_pval
MR_Pipeline_pval_cols<-colnames(decode_biogen[, grep("MR_Pipeline_pval$", colnames(decode_biogen))])

MR_Pipeline_miss<-decode_biogen[!complete.cases(decode_biogen[MR_Pipeline_pval_cols]), ]
MR_Pipeline_nomiss<-decode_biogen[complete.cases(decode_biogen[MR_Pipeline_pval_cols]), ]

MR_Pipeline_df <- MR_Pipeline_nomiss[, grep("MR_Pipeline_pval$", colnames(decode_biogen))]
metap_MR_Pipeline<-vector()

for(indx in 1:nrow(MR_Pipeline_df)) {
  teest<-unlist(as.vector(MR_Pipeline_df[indx,]))
  metap_MR_Pipeline<-c(metap_MR_Pipeline,allmetap(teest,method="sumlog")$p[[1]] )
}


MR_Pipeline_nomiss$metap_BiogenMR_Pipeline<-metap_MR_Pipeline
MR_Pipeline_miss$metap_BiogenMR_Pipeline<-"NA"
decode_biogen <- rbind(MR_Pipeline_nomiss,MR_Pipeline_miss)


##Reordering the columns
first_cols <- colnames(decode_biogen[, grep("MR_Pipeline_pval$|MRIVWtest_Pvalue$", colnames(decode_biogen))])
second_cols <- colnames(decode_biogen[, grep("_Outcome$|_Exposure$", colnames(decode_biogen))])

new_order <- c(
  c("Gene_Symbol"), first_cols, c("metap_BiogenMR_Pipeline", "metap_MRIVWtest"), second_cols,
  setdiff(colnames(decode_biogen), c("Gene_Symbol", first_cols, "metap_BiogenMR_Pipeline", "metap_MRIVWtest", second_cols))
)


decode_biogen_reordered_df <- decode_biogen[, new_order]
write.csv(decode_biogen_reordered_df, "Biogen_Decode_pQTL_trans_exposur_noMHC_CompleteMR_AnalysisResults_Meetap.csv", row.names = FALSE)

