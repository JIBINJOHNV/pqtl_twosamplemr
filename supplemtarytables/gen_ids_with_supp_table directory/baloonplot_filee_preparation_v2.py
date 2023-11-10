
import pandas as pd 
import numpy as np
import math


def to_scientific(num):
    return '{:.2e}'.format(num)

df=pd.read_excel("MR_Supplemetary_Table_New.xlsx",sheet_name="Sheet11")

df2=df[['Gene symbol', 'CIS/TRANS','Outcome','Estimate_Direction','Lowest_Pvalue','FDR_P_value','BONF_P_Value']]
df2=df2[~df2['Estimate_Direction'].isin(['+,-','-,+'])]
#df2=df2.drop("Estimate_Direction",axis=1)

replace_mapping = {'NA,+':'+',"+,+":"+",'+,NA':'+',"-,-":'-',"-,NA":'-','NA,-':'-',}
# Replace values in the 'Estimate_Direction' column
df2['Estimate_Direction'] = df2['Estimate_Direction'].replace(replace_mapping)


df2['Outcome']=df2['Outcome'].map({'BIP_PGC3_noukb':"BIP", 'Depression_iPSYCH_2023':"MDD", 'PGC3_SCZ_NoUKB':"SCZ",'Cognition_Meta':"COG"})
df2=df2.sort_values(by=["Gene symbol","CIS/TRANS","Outcome","Lowest_Pvalue"])
df2=df2.drop_duplicates(subset=["Gene symbol","CIS/TRANS","Outcome"],keep="first")

conditions = [
    ( (df2["BONF_P_Value"] < 0.05) & df2["BONF_P_Value"].notna()),
    ((df2["BONF_P_Value"] > 0.05) & (df2["FDR_P_value"] < 0.05)  & df2["FDR_P_value"].notna()  & df2["BONF_P_Value"].notna()),
    ((df2["FDR_P_value"] > 0.05) & (df2['Lowest_Pvalue'] < 0.05) & df2["FDR_P_value"].notna()  & df2["BONF_P_Value"].notna())
]

values=["strictly significant","FDR","nominally significant"]
df2["shapee_vale"]=np.select(conditions, values, default='Unknown')

df3=df2.pivot(index=['Gene symbol','CIS/TRANS'], columns='Outcome', values=['Lowest_Pvalue',"shapee_vale","BONF_P_Value","FDR_P_value","Estimate_Direction"]).reset_index()


df3.columns = ['_'.join(col).strip() for col in df3.columns.values]

fdr_columns=['FDR_P_value_BIP', 'FDR_P_value_COG', 'FDR_P_value_MDD', 'FDR_P_value_SCZ']
bnf_columns=['BONF_P_Value_BIP', 'BONF_P_Value_COG', 'BONF_P_Value_MDD', 'BONF_P_Value_SCZ']
bnf_df = df3[(df3[bnf_columns] < 0.05).any(axis=1)]



bnf_05_m1=pd.melt(bnf_df,id_vars=['Gene symbol_', 'CIS/TRANS_'],
                        value_vars=['Lowest_Pvalue_BIP', 'Lowest_Pvalue_COG', 'Lowest_Pvalue_MDD', 'Lowest_Pvalue_SCZ'])
bnf_05_m1=bnf_05_m1.rename(columns={"value":"P-value","variable":"Phenotype"})
bnf_05_m1["Phenotype"]=bnf_05_m1["Phenotype"].str.replace("Lowest_Pvalue_","")


bnf_05_m2=pd.melt(bnf_df,id_vars=['Gene symbol_', 'CIS/TRANS_'],
                        value_vars=['shapee_vale_BIP', 'shapee_vale_COG', 'shapee_vale_MDD', 'shapee_vale_SCZ'])
bnf_05_m2=bnf_05_m2.rename(columns={"value":"shape_vale","variable":"Phenotype"})
bnf_05_m2["Phenotype"]=bnf_05_m2["Phenotype"].str.replace("shapee_vale_","")


bnf_05_m3=pd.melt(bnf_df,id_vars=['Gene symbol_', 'CIS/TRANS_'],
                        value_vars=['Estimate_Direction_BIP','Estimate_Direction_COG', 'Estimate_Direction_MDD','Estimate_Direction_SCZ'])
bnf_05_m3=bnf_05_m3.rename(columns={"value":"Estimate_Direction","variable":"Phenotype"})
bnf_05_m3["Phenotype"]=bnf_05_m3["Phenotype"].str.replace("Estimate_Direction_","")


bnf_05_m=pd.merge(bnf_05_m1,bnf_05_m2,on=['Gene symbol_', 'CIS/TRANS_','Phenotype'])
bnf_05_m=pd.merge(bnf_05_m,bnf_05_m3,on=['Gene symbol_', 'CIS/TRANS_','Phenotype'],how="outer")
bnf_05_m=bnf_05_m[~bnf_05_m["shape_vale"].isna()]
bnf_05_m=bnf_05_m.rename(columns={"Gene symbol_":"Gene","CIS/TRANS_":"CIS/TRANS",'Estimate_Direction_':'Estimate_Direction'})
bnf_05_m['P-value'] = bnf_05_m['P-value'].apply(lambda x: -math.log10(x))
bnf_05_m=bnf_05_m[bnf_05_m["P-value"]>=1.300000]


bnf_05_m.to_csv("Alll_pheeotype_cis_trans_bonf.csv",index=None)
