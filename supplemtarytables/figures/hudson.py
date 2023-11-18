library(hudson)
library(dplyr)

data(phewas.t)
data(phewas.b)



cis_df=read.csv('Decode_Biogen_CisExposure_all_phenotype_hudson_input.csv')
trans_df=read.csv('Decode_Biogen_TransExposure_all_phenotype_hudson_input.csv')


negative_cis_df <- cis_df %>%filter(Estimate_Direction %in% c('-,NA', '-,-', 'NA,-'))
positive_cis_df <- cis_df %>%filter(Estimate_Direction %in% c('+,NA', '+,+', 'NA,+'))

### All phenotype but cis and trans in separate figures
negative_trans_df <- trans_df %>%filter(Estimate_Direction %in% c('-,NA', '-,-', 'NA,-'))
positive_trans_df <- trans_df %>%filter(Estimate_Direction %in% c('+,NA', '+,+', 'NA,+'))

phemirror(top=positive_cis_df, bottom =negative_cis_df, toptitle = "Beta value",opacity=1,bottomtitle = "Negative Beta value",
            annotate_p = c(0.0000135, 0.0000135),file = "cis_all_phenotype",)

phemirror(top=positive_trans_df, bottom =negative_trans_df, toptitle = "Beta value",opacity=1,bottomtitle = "Negative Beta value",
            annotate_p = c(0.00001, 0.00001),file = "trans_all_phenotype",)


combined_df <- rbind(cis_df, trans_df)
#combined_df <- subset(combined_df, select = -PHE)
#colnames(combined_df)[colnames(combined_df) == "CIS.TRANS"] <- "PHE"

## BIPOLAR disorder
bip_df=combined_df[combined_df$PHE=="BIP", ]
bip_df <- subset(bip_df, select = -PHE)
colnames(bip_df)[colnames(bip_df) == "CIS.TRANS"] <- "PHE"
negative_bip_df <- bip_df %>%filter(Estimate_Direction %in% c('-,NA', '-,-', 'NA,-'))
positive_bip_df <- bip_df %>%filter(Estimate_Direction %in% c('+,NA', '+,+', 'NA,+'))
phemirror(top=positive_bip_df, bottom =negative_bip_df,opacity=1,
            annotate_p = c(0.0000134, 0.0000134),file = "bipolar_cis_trans",toptitle = "Positive Beta value",bottomtitle = "Negative Beta value",)


#Cognition
cog_df=combined_df[combined_df$PHE=="CTP", ]
cog_df <- subset(cog_df, select = -PHE)
colnames(cog_df)[colnames(cog_df) == "CIS.TRANS"] <- "PHE"
negative_cog_df <- cog_df %>%filter(Estimate_Direction %in% c('-,NA', '-,-', 'NA,-'))
positive_cog_df <- cog_df %>%filter(Estimate_Direction %in% c('+,NA', '+,+', 'NA,+'))
phemirror(top=positive_cog_df, bottom =negative_cog_df,opacity=1,annotate_p = c(0.000010724, 0.000010724),
            file = "cognition_cis_trans",toptitle = "Positive Beta value",bottomtitle = "Negative Beta value",)

##Depression
dep_df=combined_df[combined_df$PHE=="Depression", ]
dep_df <- subset(dep_df, select = -PHE)
colnames(dep_df)[colnames(dep_df) == "CIS.TRANS"] <- "PHE"
negative_dep_df <- dep_df %>%filter(Estimate_Direction %in% c('-,NA', '-,-', 'NA,-'))
positive_dep_df <- dep_df %>%filter(Estimate_Direction %in% c('+,NA', '+,+', 'NA,+'))
phemirror(top=positive_dep_df, bottom =negative_dep_df,opacity=1,annotate_p = c(0.00001, 0.00001),
                file = "depression_cis_trans",toptitle = "Positive Beta value",bottomtitle = "Negative Beta value",)


### SCZ
scz_df=combined_df[combined_df$PHE=="SCZ", ]
scz_df <- subset(scz_df, select = -PHE)
colnames(scz_df)[colnames(scz_df) == "CIS.TRANS"] <- "PHE"

negative_scz_df <- scz_df %>%filter(Estimate_Direction %in% c('-,NA', '-,-', 'NA,-'))
positive_scz_df <- scz_df %>%filter(Estimate_Direction %in% c('+,NA', '+,+', 'NA,+'))
phemirror(top=positive_scz_df, bottom =negative_scz_df,opacity=1,annotate_p = c(0.000013, 0.000013),
                    file = "scz_cis_trans",toptitle = "Schizophrenia Positive Beta value",bottomtitle = "Schizophrenia Negative Beta value",)

