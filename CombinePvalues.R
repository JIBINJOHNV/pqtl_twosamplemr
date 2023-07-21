library(metap)

df=read.csv("SCZ_Only_Cis_exposure_CompleteMR_AnalysisResults_BiologicalPsychiatry_S5.csv")
df2<-df[complete.cases(df[c("p.value", "MRIVWtest_Pvalue")]), ]

miss<-df[!complete.cases(df[c("p.value", "MRIVWtest_Pvalue")]), ]

df3<-df2[,c("p.value","MRIVWtest_Pvalue")]


metap_values<-vector()

for(indx in 1:nrow(df3)) {
  teest<-unlist(as.vector(df3[indx,]))
  metap_values<-c(metap_values,allmetap(teest,method="sumlog")$p[[1]] )
}

df2$metap_British_BiologicalPsychiatry<-metap_values



df3<-df2[,c("p.value","MR_Pipeline_pval")]
metap_values<-vector()
for(indx in 1:nrow(df3)) {
  teest<-unlist(as.vector(df3[indx,]))
  metap_values<-c(metap_values,allmetap(teest,method="sumlog")$p[[1]] )
}

df2$metap_Biogen_BiologicalPsychiatry<-metap_values


result<-rbind.fill(df2,miss)

write.csv(result, "Testmeta_p.csv", row.names = FALSE)
