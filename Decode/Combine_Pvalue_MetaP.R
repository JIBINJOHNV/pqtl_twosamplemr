library(metap)
library(plyr)


df=read.csv("Biogen_Decode_pQTL.csv")
scolumns<-c("Decode_MRIVWtest_Pvalue","Biogen_MRIVWtest_Pvalue")
df3<-df[,scolumns]

metap_values1<-vector()

for(indx in 1:nrow(df3)) {
  teest<-unlist(as.vector(df3[indx,]))
  metap_values1<-c(metap_values1,allmetap(teest,method="sumlog")$p[[1]] )
}



scolumns<-c("Decode_MR_Pipeline_pval","Biogen_MR_Pipeline_pval")
df3<-df[,scolumns]

metap_values2<-vector()

for(indx in 1:nrow(df3)) {
  teest<-unlist(as.vector(df3[indx,]))
  metap_values2<-c(metap_values2,allmetap(teest,method="sumlog")$p[[1]] )
}

df$metap_MRIVWtest<-metap_values1
df$metap_MR_Pipeline_pval<-metap_values2


write.csv(df, "Biogen_Decode_pQTL_Meetap2.csv", row.names = FALSE)




############################################





df=read.csv("Biogen_Decode_pQTL.csv")
scolumns<-c("Decode_MRIVWtest_Pvalue","Biogen_MRIVWtest_Pvalue")

df2<-df[complete.cases(df[scolumns]), ]

miss<-df[!complete.cases(df[scolumns]), ]

df3<-df2[,scolumns]


metap_values<-vector()

for(indx in 1:nrow(df3)) {
  teest<-unlist(as.vector(df3[indx,]))
  metap_values<-c(metap_values,allmetap(teest,method="sumlog")$p[[1]] )
}

df2$metap_MRIVWtest<-metap_values


scolumns<-c("Decode_MR_Pipeline_pval","Biogen_MR_Pipeline_pval")
df3<-df2[,c("p.value","MR_Pipeline_pval")]
metap_values<-vector()
for(indx in 1:nrow(df3)) {
  teest<-unlist(as.vector(df3[indx,]))
  metap_values<-c(metap_values,allmetap(teest,method="sumlog")$p[[1]] )
}

df2$metap_MR_Pipeline_pval<-metap_values


result<-rbind.fill(df2,miss)

write.csv(result, "Biogen_Decode_pQTL_Meetap.csv", row.names = FALSE)
