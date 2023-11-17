library(hudson)
library(dplyr)



cis_df=read.csv('Decode_Biogen_CisExposure_all_phenotype_hudson_input.csv')
trans_df=read.csv('Decode_Biogen_TransExposure_all_phenotype_hudson_input.csv')

negative_cis_df <- cis_df %>%filter(Estimate_Direction %in% c('-,NA', '-,-', 'NA,-'))
positive_cis_df <- cis_df %>%filter(Estimate_Direction %in% c('+,NA', '+,+', 'NA,+'))

negative_cis_df_cog=negative_cis_df[negative_cis_df$PHE=="Cognition",]
negative_cis_df_nocog=negative_cis_df[negative_cis_df$PHE!="Cognition",]

positive_cis_df_cog=positive_cis_df[positive_cis_df$PHE=="Cognition",]
positive_cis_df_nocog=positive_cis_df[positive_cis_df$PHE!="Cognition",]

negative_cis_df<-rbind(negative_cis_df_nocog,positive_cis_df_cog)
positive_cis_df<-rbind(positive_cis_df_nocog,negative_cis_df_cog)




### All phenotype but cis and trans in separate figures
negative_trans_df <- trans_df %>%filter(Estimate_Direction %in% c('-,NA', '-,-', 'NA,-'))
positive_trans_df <- trans_df %>%filter(Estimate_Direction %in% c('+,NA', '+,+', 'NA,+'))



negative_trans_df_cog=negative_trans_df[negative_trans_df$PHE=="Cognition",]
negative_trans_df_nocog=negative_trans_df[negative_trans_df$PHE!="Cognition",]

positive_trans_df_cog=positive_trans_df[positive_trans_df$PHE=="Cognition",]
positive_trans_df_nocog=positive_trans_df[positive_trans_df$PHE!="Cognition",]

negative_trans_df<-rbind(negative_trans_df_nocog,positive_trans_df_cog)
positive_trans_df<-rbind(positive_trans_df_nocog,negative_trans_df_cog)



seed=305

phemirror(top=positive_cis_df, bottom =negative_cis_df, toptitle = "Positive Beta value",opacity=1,bottomtitle = "Negative Beta value",
            annotate_p = c(0.0000135, 0.0000135),file = "cis_all_phenotype",)

phemirror(top=positive_trans_df, bottom =negative_trans_df, toptitle = "Positive Beta value",opacity=1,bottomtitle = "Negative Beta value",
            annotate_p = c(0.00001, 0.00001),file = "trans_all_phenotype",)

combined_df <- rbind(cis_df, trans_df)

#Cognition
cog_df=combined_df[combined_df$PHE=="CTP", ]
cog_df <- subset(cog_df, select = -PHE)
colnames(cog_df)[colnames(cog_df) == "CIS.TRANS"] <- "PHE"
negative_cog_df <- cog_df %>%filter(Estimate_Direction %in% c('-,NA', '-,-', 'NA,-'))
positive_cog_df <- cog_df %>%filter(Estimate_Direction %in% c('+,NA', '+,+', 'NA,+'))
phemirror(bottom=positive_cog_df, top =negative_cog_df,opacity=1,annotate_p = c(0.0000135, 0.0000135),
            file = "cognition_cis_trans_Beta_Reversed",toptitle = "CTP Positive Beta value",bottomtitle = "CTP Negative Beta value",)
