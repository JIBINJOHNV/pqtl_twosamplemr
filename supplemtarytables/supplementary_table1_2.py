import pandas as pd 
import os,glob

def supplementary_table(Pqtl_Source,Pqtl_Source_label,pqt_type,pqt_type_label,gwas):
    basedir="/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/"
    basedir1=f"{Pqtl_Source}_CombinedResultsfromAllBatches/DiseaseSpecific/"            
    basedir2=f"{gwas}/{Pqtl_Source}_{gwas}_{pqt_type}/"
    ivdelta=f"{basedir}{basedir1}{basedir2}/{Pqtl_Source}_{gwas}_{pqt_type}_MendelianRandomization_IVW_Delta_Test.csv"
    all=f"{basedir}{basedir1}{basedir2}/{Pqtl_Source}_{gwas}_{pqt_type}_MendelianRandomization_AllTest.csv"
    direction=f"{basedir}{basedir1}{basedir2}/{Pqtl_Source}_{gwas}_{pqt_type}_TwoSampleMR_Analysis_Directionality_Test.csv"
    het=f"{basedir}{basedir1}{basedir2}/{Pqtl_Source}_{gwas}_{pqt_type}_TwoSampleMR_Analysis_HeterogenityTest.csv"
    hp=f"{basedir}{basedir1}{basedir2}/{Pqtl_Source}_{gwas}_{pqt_type}_TwoSampleMR_Analysis_Hpleiotropy_Test.csv"
    iv_df=pd.read_csv(ivdelta).rename(columns={'Exposure':"exposure",'SNPs':'N_SNPs'})[['Model', 'Outcome', 'Penalized', 'Estimate', 'CILower', 'Alpha','Heter.Stat', 'exposure', 'Robust', 'Correlation', 'StdError','CIUpper', 'Pvalue', 'RSE','N_SNPs']]
    all_df=pd.read_csv(all)[["exposure","SNPs"]].drop_duplicates()
    dir_df=pd.read_csv(direction)[['exposure','snp_r2.exposure','snp_r2.outcome', 'correct_causal_direction', 'steiger_pval']]
    hp_df=pd.read_csv(hp)[['exposure', 'egger_intercept','se', 'pval']]
    het_df=pd.read_csv(het)[['exposure', 'method', 'Q','Q_df', 'Q_pval']]
    het_df=het_df.replace({'Maximum likelihood':'Maximum_likelihood', 'MR Egger':'MR_Egger','Inverse variance weighted':"IVW",'Unweighted regression':'Unweighted_regression','IVW radial':'IVW_radial'})
    het_df=pd.pivot_table(het_df,values=["Q", "Q_df","Q_pval"],index="exposure",columns="method").reset_index()
    het_df.columns=[' '.join(col).strip() for col in het_df.columns.values]
    iv_all_df=pd.merge(iv_df,all_df,on="exposure",how="left")
    iv_all_dir_df=pd.merge(iv_all_df,dir_df,on="exposure",how="left")
    iv_all_dir_hp_df=pd.merge(iv_all_dir_df,hp_df,on="exposure",how="left")
    iv_all_dir_hp_het_df=pd.merge(iv_all_dir_hp_df,het_df,on="exposure",how="left")
    if Pqtl_Source=="Decode":
        iv_all_dir_hp_het_df["Protein"]=iv_all_dir_hp_het_df["exposure"].str.split("_",expand=True)[0]
        iv_all_dir_hp_het_df["Panel"]="NA"
    if Pqtl_Source=="Biogen":
        iv_all_dir_hp_het_df["Protein"]=iv_all_dir_hp_het_df["exposure"].str.split(":",expand=True)[0]
        iv_all_dir_hp_het_df["Panel"]=iv_all_dir_hp_het_df['exposure'].str.split(":",expand=True)[2]
    iv_all_dir_hp_het_df["CIS/TRANS"]=pqt_type_label
    iv_all_dir_hp_het_df["Pqtl_Source"]=Pqtl_Source_label
    order=['Protein','CIS/TRANS','Pqtl_Source',"Panel"]+list(iv_all_dir_hp_het_df.columns[:-4])
    iv_all_dir_hp_het_df=iv_all_dir_hp_het_df[order]
    iv_all_dir_hp_het_df.to_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/supplementary_Tablees/{gwas}_{Pqtl_Source}_{pqt_type}_suppelementary_table.csv",index=None)



### Create supplementary table 1
PqtlSource=["Decode","Biogen"]
pqtls={"CIS":"CisExposure","TRANS":"TransExposureNoMHC"}
gwases=["Depression_iPSYCH_2023" ,"PGC3_SCZ_NoUKB" ,"Cognition" , "BIP_PGC3_noukb"]
cis_trans_decode_biogen_df=pd.DataFrame()

for gwas in gwases:
    for source in PqtlSource:
        cis_file=f"/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/supplementary_Tablees/{gwas}_{source}_{pqtls['CIS']}_suppelementary_table.csv"
        trans_file=f"/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/supplementary_Tablees/{gwas}_{source}_{pqtls['TRANS']}_suppelementary_table.csv"
        cis_df=pd.read_csv(cis_file)
        cis_df=cis_df.drop(['CIS/TRANS'],axis=1)
        cis_df["CIS"]="YES"
        prefix="CIS_"
        cis_df.columns=[ prefix+x for x in cis_df.columns ]
        cis_df.rename(columns={prefix+'Protein':'Protein',prefix+'Outcome':'Outcome',prefix+'Pqtl_Source':'Pqtl_Source',prefix+'exposure':'exposure',prefix+"CIS":"CIS"},inplace=True)
        trans_df=pd.read_csv(trans_file)
        trans_df=trans_df.drop(['CIS/TRANS'],axis=1)
        trans_df["TRANS"]="YES"
        prefix="TRANS_"
        trans_df.columns=[ prefix+x for x in trans_df.columns ]
        trans_df.rename(columns={prefix+'Protein':'Protein',prefix+'Outcome':'Outcome',prefix+'Pqtl_Source':'Pqtl_Source',prefix+'exposure':'exposure',prefix+"TRANS":"TRANS"},inplace=True)
        
        cistrans_df=pd.merge(cis_df,trans_df,on=['Protein','Outcome','Pqtl_Source','exposure'],how="outer")
        cistrans_df[['TRANS', 'CIS']]=cistrans_df[['TRANS', 'CIS']].fillna("NO")
        order1=['Protein','Pqtl_Source','Outcome','TRANS', 'CIS']
        order2=[x for x in cistrans_df.columns if x not in order1 ]
        cistrans_df=cistrans_df[order1+order2]
        cis_snp_df=cistrans_df[['CIS_SNPs','Protein']]
        cis_snp_df['CIS_SNPs']=cis_snp_df['CIS_SNPs'].str.split(",")
        cis_snp_df=cis_snp_df.explode('CIS_SNPs')
        cis_snp_df=cis_snp_df[~cis_snp_df["CIS_SNPs"].isna()]
        cistrans_df['TRANS_SNPs_to_split']=cistrans_df['TRANS_SNPs'].str.split(",")
        trans_snp_df=cistrans_df[['TRANS_SNPs_to_split','TRANS_SNPs']].explode('TRANS_SNPs_to_split')
        cis_trans=pd.merge(trans_snp_df,cis_snp_df,left_on="TRANS_SNPs_to_split",right_on="CIS_SNPs")
        cis_trans=cis_trans[['TRANS_SNPs','CIS_SNPs','Protein']].drop_duplicates()
        cis_trans=cis_trans.rename(columns={'CIS_SNPs':"Implicated_CIS_SNPs_In_Trans",'Protein':"Protein_AssociatedWith_Implicated_CIS_SNPs_In_Trans"})
        cis_trans=cis_trans.groupby("TRANS_SNPs").agg({"Implicated_CIS_SNPs_In_Trans": ','.join,"Protein_AssociatedWith_Implicated_CIS_SNPs_In_Trans": ','.join}).reset_index()
        cis_trans=cis_trans.drop_duplicates()
        cistrans_df2=pd.merge(cistrans_df,cis_trans,on="TRANS_SNPs",how="left").drop("TRANS_SNPs_to_split",axis=1).drop_duplicates()
        cis_trans_decode_biogen_df=pd.concat([cis_trans_decode_biogen_df,cistrans_df2])
        cistrans_df2.to_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/supplementary_Tablees/{gwas}_{source}_CisExposure_TransExposureNoMHC_suppelementary_table.csv",index=None)


cis_trans_decode_biogen_df.to_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/supplementary_Tablees/Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table.csv",index=None)



#################   Supplementary Tabl2
PqtlSource=["Decode","Biogen"]
gwases=["Depression_iPSYCH_2023" ,"PGC3_SCZ_NoUKB" ,"Cognition" , "BIP_PGC3_noukb"]
cis_trans_decode_biogen_df=pd.DataFrame()
pqtls=["CisExposure","TransExposureNoMHC"]

        #gwas=gwases[0]          
        #pqtl=pqtls[0]   
result_df=pd.DataFrame()

for gwas in gwases:
    for pqtl in pqtls:      
        decode_cis_file=f"/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/supplementary_Tablees/{gwas}_Decode_{pqtl}_suppelementary_table.csv"
        biogeen_cis_file=f"/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/supplementary_Tablees/{gwas}_Biogen_{pqtl}_suppelementary_table.csv"
        cis_metap_file=f"/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/Metap/MetaP_AnalysisResults/Biogen_Decode_pQTL_{gwas}_{pqtl}_British_IVDelt_MetapAnalysis.csv"
        decode_cis_df=pd.read_csv(decode_cis_file)
        decode_cis_df=decode_cis_df.drop("Pqtl_Source",axis=1)
        prefix="DECODE_"
        decode_cis_df.columns=[ prefix+x for x in decode_cis_df.columns ]
        decode_cis_df.rename(columns={prefix+'Protein':'Protein',prefix+'Outcome':'Outcome',prefix+"CIS/TRANS":"CIS/TRANS"},inplace=True)
        
        biogen_cis_df=pd.read_csv(biogeen_cis_file)
        biogen_cis_df=biogen_cis_df.drop("Pqtl_Source",axis=1)
        prefix="UKBB_PPP_"
        biogen_cis_df.columns=[ prefix+x for x in biogen_cis_df.columns ]
        biogen_cis_df.rename(columns={prefix+'Protein':'Protein',prefix+'Outcome':'Outcome',prefix+"CIS/TRANS":"CIS/TRANS"},inplace=True)
        decode_biogen_cis_df=pd.merge(decode_cis_df,biogen_cis_df,on=['Protein','Outcome',"CIS/TRANS"],how="outer")
        
        cis_metap_df=pd.read_csv(cis_metap_file)[['Gene_Symbol','metap_MRIVWtest', 'outcome', 'Decode.PQTL_exposure','Biogen.PQTL_outcome']]
        cis_metap_df=cis_metap_df[~cis_metap_df["metap_MRIVWtest"].isna()]
        cis_metap_df.rename(columns={'Gene_Symbol':'Protein', 'metap_MRIVWtest':'metap_Pvalue', 'outcome':'Outcome', 
                                    'Decode.PQTL_exposure':'DECODE_exposure','Biogen.PQTL_outcome':'UKBB_PPP_exposure'},inplace=True)
        dbm_cis_df=pd.merge(decode_biogen_cis_df,cis_metap_df,on=['Protein','Outcome','DECODE_exposure','UKBB_PPP_exposure'],how="left")
        conditions = [((~dbm_cis_df['DECODE_exposure'].isna()) & (~dbm_cis_df['UKBB_PPP_exposure'].isna())),
                    (dbm_cis_df['DECODE_exposure'].isna()) & (~dbm_cis_df['UKBB_PPP_exposure'].isna()),
                    (~dbm_cis_df['DECODE_exposure'].isna()) & (dbm_cis_df['UKBB_PPP_exposure'].isna())]
        
        values = ['BOTH', 'UKBB-PPP_Specific', 'DECODE_Specific']
        dbm_cis_df['Cateegary'] = np.select(conditions, values)
        dbm_cis_df["Lowest_Pvalue"]=dbm_cis_df[['DECODE_Pvalue', 'UKBB_PPP_Pvalue', 'metap_Pvalue']].apply(lambda row: row.min(), axis=1)
        dbm_cis_df["Lowest_Pvalue_Source"]=dbm_cis_df[['DECODE_Pvalue', 'UKBB_PPP_Pvalue', 'metap_Pvalue']].idxmin(axis=1).str.replace("_Pvalue","")
        dbm_cis_df["Largest_Pvalue"]=dbm_cis_df[['DECODE_Pvalue', 'UKBB_PPP_Pvalue', 'metap_Pvalue']].apply(lambda row: row.max(), axis=1)
        dbm_cis_df["Largest_Pvalue_Source"]=dbm_cis_df[['DECODE_Pvalue', 'UKBB_PPP_Pvalue', 'metap_Pvalue']].idxmax(axis=1).str.replace("_Pvalue","")
        
        conditions2 = [(dbm_cis_df['DECODE_Estimate'].isna() & (dbm_cis_df['UKBB_PPP_Estimate'] >= 0)).astype(bool),
                ((dbm_cis_df['DECODE_Estimate'] >= 0) & dbm_cis_df['UKBB_PPP_Estimate'].isna()).astype(bool),
                (dbm_cis_df['DECODE_Estimate'].isna() & (dbm_cis_df['UKBB_PPP_Estimate'] < 0)).astype(bool),
                ((dbm_cis_df['DECODE_Estimate'] < 0) & dbm_cis_df['UKBB_PPP_Estimate'].isna()).astype(bool),
                ((dbm_cis_df['DECODE_Estimate'] < 0) & (dbm_cis_df['UKBB_PPP_Estimate'] < 0)).astype(bool),
                ((dbm_cis_df['DECODE_Estimate'] >= 0) & (dbm_cis_df['UKBB_PPP_Estimate'] >= 0)).astype(bool),
                ((dbm_cis_df['DECODE_Estimate'] >= 0) & (dbm_cis_df['UKBB_PPP_Estimate'] < 0)).astype(bool),
                ((dbm_cis_df['DECODE_Estimate'] < 0) & (dbm_cis_df['UKBB_PPP_Estimate'] >= 0)).astype(bool)]
        
        values2 = ["NA,+", "+,NA", "NA,-", "-,NA", "-,-", "+,+", "+,-", "-,+"]
        dbm_cis_df['Estimate_Direction'] = np.select(conditions2, values2, default="Default_Value")
        order1=['Protein', 'CIS/TRANS','Outcome','Cateegary','Estimate_Direction','metap_Pvalue','Lowest_Pvalue', 'Largest_Pvalue','Lowest_Pvalue_Source','Largest_Pvalue_Source']
        order2=['DECODE_Panel','DECODE_exposure','UKBB_PPP_exposure','UKBB_PPP_Panel']
        order3=[x for x in dbm_cis_df.columns if x not in order1 and x not in order2]
        order=order1+order3+order2
        dbm_cis_df=dbm_cis_df[order]
        result_df=pd.concat([result_df,dbm_cis_df])


result_df.to_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/supplementary_Tablees/Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2.csv",index=None)

df=pd.read_csv("/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/supplementary_Tablees/Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table.csv")
df2=df[['Protein','TRANS_SNPs','Implicated_CIS_SNPs_In_Trans','Protein_AssociatedWith_Implicated_CIS_SNPs_In_Trans']]

result_df2=result_df[['Protein','CIS/TRANS','DECODE_SNPs','UKBB_PPP_SNPs']]

merge1=pd.merge(result_df2,df2,left_on=["Protein","DECODE_SNPs"],right_on=['Protein', 'TRANS_SNPs'])
merge1=merge1[~merge1["Implicated_CIS_SNPs_In_Trans"].isna()]
merge2=pd.merge(result_df2,df2,left_on=["Protein","UKBB_PPP_SNPs"],right_on=['Protein', 'TRANS_SNPs'])
merge2=merge2[~merge2["Implicated_CIS_SNPs_In_Trans"].isna()]

merge=pd.concat([merge1,merge2]).drop('TRANS_SNPs',axis=1).drop_duplicates()

result_df5=pd.merge(result_df,merge,on=['Protein', 'CIS/TRANS',"DECODE_SNPs","UKBB_PPP_SNPs"],how="left")
result_df5.to_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/supplementary_Tablees/Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2b.csv",index=None)




