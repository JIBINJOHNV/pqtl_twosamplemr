import pandas as pd

df=pd.read_csv("Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2_withGeneInfo.csv")
df2=df[['Approved symbol', 'CIS/TRANS','Outcome','Estimate_Direction','Lowest_Pvalue']]
df2=df2[~df2['Estimate_Direction'].isin(['+,-','-,+'])]

df2=df2.drop("Estimate_Direction",axis=1)


df2['Outcome']=df2['Outcome'].map({'BIP_PGC3_noukb':"BIP", 'Depression_iPSYCH_2023':"Depression", 'PGC3_SCZ_NoUKB':"SCZ",'Cognition_Meta':"Cognition"})

df2=df2.sort_values(by=["Approved symbol","CIS/TRANS","Outcome","Lowest_Pvalue"])
df2=df2.drop_duplicates(subset=["Approved symbol","CIS/TRANS","Outcome"],keep="first")
df3=df2.pivot(index=['Approved symbol','CIS/TRANS'], columns='Outcome', values='Lowest_Pvalue').reset_index()

selected_rows_0005 = df3[(df3[['BIP', 'Cognition', 'Depression', 'SCZ']] < 0.0005).any(axis=1)]
selected_rows_00005 = df3[(df3[['BIP', 'Cognition', 'Depression', 'SCZ']] < 0.00005).any(axis=1)]
selected_rows_000005 = df3[(df3[['BIP', 'Cognition', 'Depression', 'SCZ']] < 0.000005).any(axis=1)]


# Assuming 'balloon2' is your DataFrame
columns_to_check = ['BIP', 'Cognition', 'Depression', 'SCZ']
threshold = 0.01

# Replace values less than 0.01 with NaN
selected_rows_000005[columns_to_check] = selected_rows_000005[columns_to_check].mask(selected_rows_000005[columns_to_check] >threshold)
selected_rows_000005=pd.melt(selected_rows_000005, id_vars=['Approved symbol', 'CIS/TRANS'], value_vars=columns_to_check, var_name='p value')

selected_rows_000005=selected_rows_000005.rename(columns={"Approved symbol":"Grne","CIS/TRANS":"CIS/TRANS","p value":"Phenotype","value":"P-value"})
selected_rows_000005["P-value"]=-np.log10(selected_rows_000005["P-value"])
selected_rows_000005.to_csv("Alll_pheeotype_cis_trans_000005.csv",index=None)
