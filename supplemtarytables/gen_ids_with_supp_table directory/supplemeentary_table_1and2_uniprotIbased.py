
import pandas as pd 
import os,glob

### Create supplementary table 1
PqtlSource=["Decode","Biogen"]
pqtls={"CIS":"CisExposure","TRANS":"TransExposureNoMHC"}
gwases=["Depression_iPSYCH_2023" ,"PGC3_SCZ_NoUKB" ,"Cognition" , "BIP_PGC3_noukb"]
cis_trans_decode_biogen_df=pd.DataFrame()

meta_df=pd.read_csv("/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/Complete_UKBB-PPP_Deecode_protein_withcomplete_Geneinformation_SeelecteedColumns_Annotation_Drug_cellstr_size_striped.csv")
meta_df[['HGNC', 'NCBI_Gene_ID','Gene_start', 'Gene_end']]=meta_df[['HGNC', 'NCBI_Gene_ID','Gene_start', 'Gene_end']].fillna(0).astype("int")
meta_df['UniProt_ID']=meta_df['UniProt_ID'].astype("str").str.strip()
meta_df=meta_df[meta_df["CHR"]!='HSCHR19_4_CTG3_1']
meta_df=meta_df[meta_df["CHR"]!='HSCHR17_7_CTG4']

update_data = [
    {"UniProt_ID": "Q14160", "CHR": 8, 'Gene_start': 143790920, 'Gene_end': 143815773},
    {"UniProt_ID": "Q8IXS6", "CHR": 9, 'Gene_start': 109640788, 'Gene_end': 109946703},
    {"UniProt_ID": "O15389", "CHR": 19, 'Gene_start': 51630101, 'Gene_end': 51645545},
    {"UniProt_ID": "O43423", "CHR": 4, 'Gene_start': 164197007, 'Gene_end': 164197711}]

# Create a dictionary to map UniProt_ID to the update data
update_dict = {item['UniProt_ID']: item for item in update_data}
# Update the DataFrame row by row
for index, row in meta_df.iterrows():
    uni_id = row['UniProt_ID']
    if uni_id in update_dict:
        update_row = update_dict[uni_id]
        meta_df.at[index, 'CHR'] = update_row['CHR']
        meta_df.at[index, 'Gene_start'] = update_row['Gene_start']
        meta_df.at[index, 'Gene_end'] = update_row['Gene_end']



# Your dictionary
C2orf66 = {'DECODE_Protein':'C2orf66','DECODE_exposure':'C2orf66_5677_15',
    'DECODE_seqid':'5677-15','DECODE_platform':'SomaScan v4',
    'DECODE_platform.1':'SomaScan v4','DECODE_target_full_name':'Chromosome 2 Open Reading Frame 66',
    'DECODE_target_name':'C2orf66','DECODE_gene_name':'C2orf66',
    'DECODE_uniprot':'Q6UXQ4','CHR':2,'Gene_start':196804417,'Gene_end':196809355,
    'Ensembl_gene_ID':'ENSG00000187944', 'NCBI_Gene_ID':401027,'HGNC':33809,'UniProt_ID':'Q6UXQ4', 'Approved symbol':'C2orf66'}

# Create a DataFrame from the dictionary
df = pd.DataFrame([C2orf66])

meta_df=pd.concat([meta_df,df])


decode_mdf=meta_df[['DECODE_exposure','UniProt_ID']].drop_duplicates()
decode_mdf=decode_mdf[~decode_mdf['DECODE_exposure'].isna()]
decode_mdf=decode_mdf.drop_duplicates()
biogen_mdf=meta_df[['UKBB_PPP_exposure','UniProt_ID']].drop_duplicates()
biogen_mdf=biogen_mdf[~biogen_mdf['UKBB_PPP_exposure'].isna()]
biogen_mdf=biogen_mdf.drop_duplicates()

def supplementary_table(Pqtl_Source,pqt_type,pqt_type_label,gwas,decode_mdf,biogen_mdf):
    basedir="/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/"
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
    iv_all_dir_hp_het_df["CIS/TRANS"]=pqt_type_label
    iv_all_dir_hp_het_df["Pqtl_Source"]=Pqtl_Source
    if Pqtl_Source=="Decode":
        iv_all_dir_hp_het_df["Protein"]=iv_all_dir_hp_het_df["exposure"].str.split("_",expand=True)[0]
    if Pqtl_Source=="Biogen":
        iv_all_dir_hp_het_df["Protein"]=iv_all_dir_hp_het_df["exposure"].str.split(":",expand=True)[0]
    order=['Protein','CIS/TRANS','Pqtl_Source']+list(iv_all_dir_hp_het_df.columns[:-3])
    iv_all_dir_hp_het_df=iv_all_dir_hp_het_df[order]
    if Pqtl_Source=="Decode":
        iv_all_dir_hp_het_df=pd.merge(iv_all_dir_hp_het_df,decode_mdf,left_on="exposure",right_on="DECODE_exposure",how="left").drop_duplicates()
    if Pqtl_Source=="Biogen":
        iv_all_dir_hp_het_df=pd.merge(iv_all_dir_hp_het_df,biogen_mdf,left_on="exposure",right_on="UKBB_PPP_exposure",how="left").drop_duplicates()
    iv_all_dir_hp_het_df.to_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/supplementary_Tablees/{gwas}_{Pqtl_Source}_{pqt_type}_suppelementary_table.csv",index=None)



PqtlSource=["Decode","Biogen"]
pqtls={"CIS":"CisExposure","TRANS":"TransExposureNoMHC"}
gwases=["Depression_iPSYCH_2023" ,"PGC3_SCZ_NoUKB" ,"Cognition" , "BIP_PGC3_noukb"]
cis_trans_decode_biogen_df=pd.DataFrame()

for Pqtl_Source in PqtlSource:
    for gwas in gwases:
        for pqtl in pqtls.items():
            pqt_type=pqtl[1]
            pqt_type_label=pqtl[0]
            supplementary_table(Pqtl_Source,pqt_type,pqt_type_label,gwas,decode_mdf,biogen_mdf)

for gwas in gwases:
    for source in PqtlSource:
        cis_file=f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/supplementary_Tablees/{gwas}_{source}_{pqtls['CIS']}_suppelementary_table.csv"
        trans_file=f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/supplementary_Tablees/{gwas}_{source}_{pqtls['TRANS']}_suppelementary_table.csv"
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
        cistrans_df2.to_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/supplementary_Tablees/{gwas}_{source}_CisExposure_TransExposureNoMHC_suppelementary_table.csv",index=None)


cis_trans_decode_biogen_df["UniProt_ID"]=np.where(cis_trans_decode_biogen_df["TRANS_UniProt_ID"].isna(),cis_trans_decode_biogen_df['CIS_UniProt_ID'],cis_trans_decode_biogen_df["TRANS_UniProt_ID"])
cis_trans_decode_biogen_df=cis_trans_decode_biogen_df.drop(["TRANS_UniProt_ID",'CIS_UniProt_ID'],axis=1)
cis_trans_decode_biogen_df.to_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/supplementary_Tablees/Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table.csv",index=None)

meta_df2=meta_df.drop(['DECODE_Protein', 'DECODE_exposure', 'DECODE_seqid', 'DECODE_platform','DECODE_platform.1', 'DECODE_target_full_name', 'DECODE_target_name','DECODE_gene_name', 'DECODE_uniprot','UKBB_PPP_Protein', 'UKBB_PPP_Panel', 'UKBB_PPP_exposure','UKBB_PPP_oid', 'UKBB_PPP_uniprot', 'UKBB_PPP_target_full_name','UKBB_PPP_gene_name', 'UKBB_PPP_platform', 'UKBB_PPP_panel'],axis=1).drop_duplicates()
result1=pd.merge(cis_trans_decode_biogen_df,meta_df2,on="UniProt_ID",how="left").drop_duplicates()
result1.to_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/supplementary_Tablees/Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table_withGeneInfo.csv",index=None)


######### Suppeleeentary table2

#################   Supplementary Tabl2
PqtlSource=["Decode","Biogen"]
gwases=["Depression_iPSYCH_2023" ,"PGC3_SCZ_NoUKB" ,"Cognition" , "BIP_PGC3_noukb"]
cis_trans_decode_biogen_df=pd.DataFrame()
pqtls=["CisExposure","TransExposureNoMHC"]

pqtls={"CIS":"CisExposure","TransExposureNoMHC":"TransExposureNoMHC"}

        #gwas=gwases[0]          
        #pqtl=pqtls[0]   

result_df=pd.DataFrame()


for gwas in gwases:
    for pre_pqtl in pqtls.items():
        pqtl=pre_pqtl[1]
        pqtl2=pre_pqtl[0]
        
        decode_cis_file=f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/supplementary_Tablees/{gwas}_Decode_{pqtl}_suppelementary_table.csv"
        biogeen_cis_file=f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/supplementary_Tablees/{gwas}_Biogen_{pqtl}_suppelementary_table.csv"
        cis_metap_file=f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/ForMeta/MetaP_AnalysisResults/Biogen_Decode_pQTL_{gwas}_{pqtl}_British_IVDelt_MetapAnalysis.csv"
        
        decode_cis_df=pd.read_csv(decode_cis_file)
        decode_cis_df=decode_cis_df.drop("Pqtl_Source",axis=1)
        prefix="DECODE_"
        decode_cis_df.columns=[ prefix+x for x in decode_cis_df.columns ]
        decode_cis_df=decode_cis_df.drop('DECODE_DECODE_exposure',axis=1)
        decode_cis_df.rename(columns={prefix+'Outcome':'Outcome',prefix+"CIS/TRANS":"CIS/TRANS",'DECODE_UniProt_ID':'UniProt_ID'},inplace=True)
        decode_cis_df=decode_cis_df.drop_duplicates()
        meta_df2=meta_df[['UniProt_ID']+[x for x in meta_df.columns if "DECODE_" in x ]].drop_duplicates()
        meta_df2['DECODE_exposure']=meta_df2['DECODE_exposure'].fillna("NA")
        meta_df2 = meta_df2[meta_df2['DECODE_exposure']!="NA"].drop_duplicates()
        meta_df2.drop("DECODE_Protein",axis=1,inplace=True)
        decode_cis_df=pd.merge(decode_cis_df,meta_df2,left_on=['DECODE_exposure','UniProt_ID'],right_on=['DECODE_exposure','UniProt_ID'],how="left").drop_duplicates()
        
        biogen_cis_df=pd.read_csv(biogeen_cis_file)
        biogen_cis_df=biogen_cis_df.drop("Pqtl_Source",axis=1)
        prefix="UKBB_PPP_"
        biogen_cis_df.columns=[ prefix+x for x in biogen_cis_df.columns ]
        biogen_cis_df=biogen_cis_df.drop('UKBB_PPP_UKBB_PPP_exposure',axis=1)
        biogen_cis_df=biogen_cis_df.drop_duplicates()
        biogen_cis_df.rename(columns={prefix+'Outcome':'Outcome',prefix+"CIS/TRANS":"CIS/TRANS",'UKBB_PPP_UniProt_ID':'UniProt_ID','UKBB_PPP_UKBB_PPP_exposure':'UKBB_PPP_exposure'},inplace=True)
        meta_df3=meta_df[['UniProt_ID']+[x for x in meta_df.columns if "UKBB_PPP_" in x]].drop_duplicates()
        meta_df3['UKBB_PPP_exposure']=meta_df3['UKBB_PPP_exposure'].fillna("NA")
        meta_df3 = meta_df3[meta_df3['UKBB_PPP_exposure']!="NA"].drop_duplicates()
        meta_df3.drop("UKBB_PPP_Protein",axis=1,inplace=True)
        biogen_cis_df=pd.merge(biogen_cis_df,meta_df3,left_on=['UKBB_PPP_exposure','UniProt_ID'],right_on=['UKBB_PPP_exposure','UniProt_ID']).drop_duplicates()
        decode_biogen_cis_df=pd.merge(decode_cis_df,biogen_cis_df,on=['UniProt_ID','Outcome',"CIS/TRANS"],how="outer")
        decode_biogen_cis_df=decode_biogen_cis_df.drop_duplicates()
        
        cis_metap_df=pd.read_csv(cis_metap_file)
        cis_metap_df=cis_metap_df[['UniProt_ID', 'outcome',f'{pqtl2}_Decode.PQTL_exposure',f'{pqtl2}_Biogen.PQTL_exposure',f'{pqtl2}_{gwas}_metap_MRIVWtest']]
        cis_metap_df.columns=[x.replace(f'{pqtl2}_', '') for x in cis_metap_df.columns]
        cis_metap_df.columns=[x.replace(f'{gwas}_', '') for x in cis_metap_df.columns]
        cis_metap_df=cis_metap_df[~cis_metap_df["metap_MRIVWtest"].isna()]
        cis_metap_df.rename(columns={'metap_MRIVWtest':'metap_Pvalue', 'outcome':'Outcome','Decode.PQTL_exposure':'DECODE_exposure','Biogen.PQTL_exposure':'UKBB_PPP_exposure'},inplace=True)
        decode_biogen_cis_df[['UniProt_ID','Outcome','DECODE_exposure','UKBB_PPP_exposure']]=decode_biogen_cis_df[['UniProt_ID','Outcome','DECODE_exposure','UKBB_PPP_exposure']].fillna("NA")
        cis_metap_df[['UniProt_ID','Outcome','DECODE_exposure','UKBB_PPP_exposure']]=cis_metap_df[['UniProt_ID','Outcome','DECODE_exposure','UKBB_PPP_exposure']].fillna("NA")
        dbm_cis_df=pd.merge(decode_biogen_cis_df,cis_metap_df,on=['UniProt_ID','Outcome','DECODE_exposure','UKBB_PPP_exposure'],how="left").drop_duplicates()        
        
        conditions = [((~dbm_cis_df['DECODE_Estimate'].isna()) & (~dbm_cis_df['UKBB_PPP_Estimate'].isna())),
                            (dbm_cis_df['DECODE_Estimate'].isna()) & (~dbm_cis_df['UKBB_PPP_Estimate'].isna()),
                            (~dbm_cis_df['DECODE_Estimate'].isna()) & (dbm_cis_df['UKBB_PPP_Estimate'].isna())]
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
        
        dbm_cis_df[['UniProt_ID','Outcome','DECODE_exposure','UKBB_PPP_exposure']]=dbm_cis_df[['UniProt_ID','Outcome','DECODE_exposure','UKBB_PPP_exposure']].fillna("NA")
        meta_df5=meta_df[[x for x in meta_df.columns if "DECODE_" not in x and "UKBB_PPP_" not in x]]
        meta_df5[['UniProt_ID']]=meta_df5[['UniProt_ID']].fillna("NA")
        meta_df5=meta_df5.drop_duplicates()
        
        dbm_cis_df=pd.merge(dbm_cis_df,meta_df5,on=['UniProt_ID'],how="left").drop_duplicates()
        order1=['Approved symbol','UniProt_ID','HGNC', 'CIS/TRANS','Outcome','Cateegary','Estimate_Direction','metap_Pvalue','Lowest_Pvalue', 'Largest_Pvalue','Lowest_Pvalue_Source','Largest_Pvalue_Source']
        order2=['DECODE_exposure','UKBB_PPP_exposure']
        order3=[x for x in dbm_cis_df.columns if x not in order1 and x not in order2]
        order=order1+order3+order2
        dbm_cis_df=dbm_cis_df[order]
        result_df=pd.concat([result_df,dbm_cis_df])


result_df['Approved symbol']=np.where(result_df['Approved symbol'].isna(),result_df['UniProt_ID'],result_df['Approved symbol'])
result_df.to_csv(f"/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/supplementary_Tablees/Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2_withGeneInfo.csv",index=None)



###########---------------------------------------------------------HUDSON PLOT ---------------------------------------------------------------###############

result_2df=result_df[['Approved symbol','Outcome','CIS/TRANS','Estimate_Direction','Cateegary','Lowest_Pvalue','DECODE_Estimate','UKBB_PPP_Estimate','CHR', 'Gene_start', 'Gene_end']].drop_duplicates()
result_2df["status"]=np.where(result_2df['Lowest_Pvalue']<0.00005,'significant',"not_significant")

update_data2=[{'Approved symbol':'LGALS7_LGALS7B',"CHR":19,"Gene_start":38770968,"Gene_end":38773517}, 
 {'Approved symbol':'AMY1A_AMY1B_AMY1C',"CHR":1,"Gene_start":103655519,"Gene_end":103664554}, 
 {'Approved symbol':'BOLA2_BOLA2B',"CHR":16,"Gene_start":29453588,"Gene_end":29454964},
 {'Approved symbol':'CGB3_CGB5_CGB8',"CHR":19,"Gene_start":49022869,"Gene_end":49024333}, 
 {'Approved symbol':'NTproBNP',"CHR":1,"Gene_start":11857464,"Gene_end":11858945 }, 
 {'Approved symbol':'DEFB4A_DEFB4B',"CHR":8,"Gene_start":7894677,"Gene_end":7896716}, 
 {'Approved symbol':'EBI3_IL27',"CHR":19,"Gene_start":4229523,"Gene_end":4237528},
 {'Approved symbol':'FUT3_FUT5',"CHR":19,"Gene_start":5842888,"Gene_end":5857122},
 {'Approved symbol':'MICB_MICA',"CHR":6,"Gene_start":31494918,"Gene_end":31511124},
 {'Approved symbol':'DEFA1_DEFA1B',"CHR":8,"Gene_start":6977649,"Gene_end":6980092},
 {'Approved symbol':'DEFB104A_DEFB104B',"CHR":8,"Gene_start":7836436,"Gene_end":7841242},
 {'Approved symbol':'CKMT1A_CKMT1B',"CHR":15,"Gene_start":43692786,"Gene_end":43699222},
 {'Approved symbol':'DEFB103A_DEFB103B', "CHR":8,"Gene_start":7881392,"Gene_end":7882663},
 {'Approved symbol':'IL12A_IL12B',"CHR":3,"Gene_start":159988835,"Gene_end":159996019}]

update_data2_df=pd.DataFrame(update_data2)


result_2df_na=result_2df[result_2df['CHR'].isna()]
result_2df_notna=result_2df[~result_2df['CHR'].isna()]
result_2df_na=result_2df_na.drop(["CHR","Gene_start","Gene_end"],axis=1)

result_2df_na=pd.merge(update_data2_df,result_2df_na,on="Approved symbol")
result_2df=pd.concat([result_2df_notna,result_2df_na,])

result_2df.to_csv("Decode_Biogen_CisExposure_TransExposureNoMHC_for_ploting.csv",index=None)


result_2df['Outcome']=result_2df['Outcome'].replace({'Depression_iPSYCH_2023':'Depression', 'PGC3_SCZ_NoUKB':'SCZ', 'Cognition_Meta':'Cognition','BIP_PGC3_noukb':'BIP'})

huston=result_2df[['Outcome','Approved symbol','CHR','Gene_start','Lowest_Pvalue','CIS/TRANS','Cateegary','Estimate_Direction','DECODE_Estimate','UKBB_PPP_Estimate']].drop_duplicates()
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
