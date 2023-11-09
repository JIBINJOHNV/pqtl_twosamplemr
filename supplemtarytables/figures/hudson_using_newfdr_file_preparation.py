
import pandas as pd
import numpy as np


meta_df=pd.read_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/supplementary_Tablees/Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table_withGeneInfo.csv")
meta_df=meta_df[['UniProt_ID', 'Approved symbol','CHR', 'Gene_start']].drop_duplicates()

result_2df=pd.read_excel("MR_Supplemetary_Table_New.xlsx",sheet_name="Sheet11")
result_2df=result_2df[['Outcome',"Gene symbol",'UniProt_ID','Lowest_Pvalue','CIS/TRANS','Cateegary','Estimate_Direction','DECODE_Estimate','UKBB_PPP_Estimate']].drop_duplicates()
result_2df=pd.merge(result_2df,meta_df,on="UniProt_ID").drop_duplicates()

result_2df['Outcome']=result_2df['Outcome'].replace({'Depression_iPSYCH_2023':'Depression', 'PGC3_SCZ_NoUKB':'SCZ', 'Cognition_Meta':'Cognition','BIP_PGC3_noukb':'BIP'})
huston=result_2df[['Outcome',"Gene symbol",'CHR','Gene_start','Lowest_Pvalue','CIS/TRANS','Cateegary','Estimate_Direction','DECODE_Estimate','UKBB_PPP_Estimate']].drop_duplicates()
huston=huston.rename(columns={"Gene symbol":'Approved symbol'})


huston=huston[ (huston['Estimate_Direction']!="-,+") & (huston['Estimate_Direction']!="+,-") ]
huston=huston.sort_values(by=["Outcome","Approved symbol","CIS/TRANS","Lowest_Pvalue"])
huston=huston.drop_duplicates(subset=['Outcome','Approved symbol','CIS/TRANS'],keep="first")

huston=huston.rename(columns={'Outcome':'PHE','Approved symbol':'SNP','CHR':'CHR','Gene_start':'POS','Lowest_Pvalue':'pvalue','Cateegary':"Shape"})
huston_cis=huston[huston['CIS/TRANS']=="CIS"]
huston_trans=huston[huston['CIS/TRANS']!="CIS"]

huston_cis.to_csv('Decode_Biogen_CisExposure_all_phenotype_hudson_input.csv',index=None)
huston_trans.to_csv('Decode_Biogen_TransExposure_all_phenotype_hudson_input.csv',index=None)



###
result_3df=result_2df[ (result_2df['Estimate_Direction']!="-,+") & (result_2df['Estimate_Direction']!="+,-") ]
result_3df_neg=result_3df[result_3df['Estimate_Direction'].isin(['-,NA', '-,-','NA,-']) ].sort_values(by=["Approved symbol","Outcome","CIS/TRANS","Lowest_Pvalue"])
result_3df_pos=result_3df[result_3df['Estimate_Direction'].isin(['+,NA', '+,+','NA,+']) ].sort_values(by=["Approved symbol","Outcome","CIS/TRANS","Lowest_Pvalue"])



result_3df_neg.sort_values(by=["Approved symbol","Outcome","CIS/TRANS","Lowest_Pvalue"])

huston=result_3df_neg[['Outcome','Approved symbol','CHR','Gene_start','Lowest_Pvalue','CIS/TRANS','Estimate_Direction']].drop_duplicates()
huston=huston.rename(columns={'Outcome':'PHE','Approved symbol':'SNP','CHR':'CHR','Gene_start':'POS','Lowest_Pvalue':'pvalue'})
huston_cis=huston[huston['CIS/TRANS']=="CIS"]
huston_trans=huston[huston['CIS/TRANS']!="CIS"]
huston_cis.to_csv('Decode_Biogen_CisExposure_all_phenotype_hudson_input_negative.csv',index=None)
huston_trans.to_csv('Decode_Biogen_TransExposure_all_phenotype_hudson_input_negative.csv',index=None)



huston=result_3df_pos[['Outcome','Approved symbol','CHR','Gene_start','Lowest_Pvalue','CIS/TRANS','Estimate_Direction']].drop_duplicates()
huston=huston.rename(columns={'Outcome':'PHE','Approved symbol':'SNP','CHR':'CHR','Gene_start':'POS','Lowest_Pvalue':'pvalue'})
huston_cis=huston[huston['CIS/TRANS']=="CIS"]
huston_trans=huston[huston['CIS/TRANS']!="CIS"]
huston_cis.to_csv('Decode_Biogen_CisExposure_all_phenotype_hudson_input_positive.csv',index=None)
huston_trans.to_csv('Decode_Biogen_TransExposure_all_phenotype_hudson_input_positive.csv',index=None)
