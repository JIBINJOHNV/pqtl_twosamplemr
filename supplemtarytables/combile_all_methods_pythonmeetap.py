

import pandas as pd 
from scipy.stats import combine_pvalues
import numpy as np

uk_common_columns=['Approved symbol','UKBB_PPP_Protein','UniProt_ID','HGNC','UKBB_PPP_exposure','Gene_Symbol','outcome','exposure','CIS/TRANS']

TwoSampleMR=['TwoSampleMR_nsnp_Ivw','TwoSampleMR_nsnp_WaldRatio','BritishMR_SNPs','BritishMR_nSNPs',
            'TwoSampleMR_Pvalue_Ivw','TwoSampleMR_Beta_Ivw', 'TwoSampleMR_SE_Ivw', 
            'TwoSampleMR_Pvalue_Ivw(FixedEffects)','TwoSampleMR_Beta_Ivw(FixedEffects)','TwoSampleMR_SE_Ivw(FixedEffects)',
            'TwoSampleMR_Pvalue_Ivw(M-Randomeffects)','TwoSampleMR_Beta_Ivw(M-Randomeffects)','TwoSampleMR_SE_Ivw(M-Randomeffects)',
            'TwoSampleMR_Pvalue_IvwRadial','TwoSampleMR_Beta_IvwRadial','TwoSampleMR_SE_IvwRadial',
            'TwoSampleMR_Pvalue_MaximumLikelihood', 'TwoSampleMR_Beta_MaximumLikelihood','TwoSampleMR_SE_MaximumLikelihood',
            'TwoSampleMR_Pvalue_MrEgger','TwoSampleMR_Beta_MrEgger','TwoSampleMR_SE_MrEgger',
            'TwoSampleMR_Pvalue_MrEgger(Bootstrap)', 'TwoSampleMR_Beta_MrEgger(Bootstrap)','TwoSampleMR_SE_MrEgger(Bootstrap)',
            'TwoSampleMR_Pvalue_PenalisedWeightedMedian', 'TwoSampleMR_Beta_PenalisedWeightedMedian', 'TwoSampleMR_SE_PenalisedWeightedMedian',
            'TwoSampleMR_Pvalue_SignConcordanceTest','TwoSampleMR_Beta_SignConcordanceTest','TwoSampleMR_SE_SignConcordanceTest',
            'TwoSampleMR_Pvalue_SimpleMedian','TwoSampleMR_Beta_SimpleMedian', 'TwoSampleMR_SE_SimpleMedian',
            'TwoSampleMR_Pvalue_SimpleMode', 'TwoSampleMR_Beta_SimpleMode','TwoSampleMR_SE_SimpleMode',
            'TwoSampleMR_Pvalue_SimpleMode(Nome)','TwoSampleMR_Beta_SimpleMode(Nome)','TwoSampleMR_SE_SimpleMode(Nome)',
            'TwoSampleMR_Pvalue_UnweightedRegression','TwoSampleMR_Beta_UnweightedRegression','TwoSampleMR_SE_UnweightedRegression',
            'TwoSampleMR_Pvalue_WaldRatio','TwoSampleMR_Beta_WaldRatio','TwoSampleMR_SE_WaldRatio', 
            'TwoSampleMR_Pvalue_WeightedMedian','TwoSampleMR_Beta_WeightedMedian','TwoSampleMR_SE_WeightedMedian', 
            'TwoSampleMR_Pvalue_WeightedMode','TwoSampleMR_Beta_WeightedMode','TwoSampleMR_SE_WeightedMode',
            'TwoSampleMR_Pvalue_WeightedMode(Nome)','TwoSampleMR_Beta_WeightedMode(Nome)','TwoSampleMR_SE_WeightedMode(Nome)']

BritishMR=['BritishMR_Pvalue_Ivw','BritishMR_Beta_Ivw','BritishMR_SE_Ivw',
            'BritishMR_Pvalue_Mr-Egger','BritishMR_Beta_Mr-Egger', 'BritishMR_SE_Mr-Egger', 
            'BritishMR_Pvalue_PenalizedIvw','BritishMR_Beta_PenalizedIvw', 'BritishMR_SE_PenalizedIvw', 
            'BritishMR_Pvalue_PenalizedMr-Egger','BritishMR_Beta_PenalizedMr-Egger','BritishMR_SE_PenalizedMr-Egger', 
            'BritishMR_Pvalue_PenalizedRobustIvw','BritishMR_Beta_PenalizedRobustIvw','BritishMR_SE_PenalizedRobustIvw',
            'BritishMR_Pvalue_PenalizedRobustMr-Egger','BritishMR_Beta_PenalizedRobustMr-Egger','BritishMR_SE_PenalizedRobustMr-Egger', 
            'BritishMR_Pvalue_PenalizedWeightedMedian','BritishMR_Beta_PenalizedWeightedMedian','BritishMR_SE_PenalizedWeightedMedian',
            'BritishMR_Pvalue_RobustIvw','BritishMR_Beta_RobustIvw','BritishMR_SE_RobustIvw',
            'BritishMR_Pvalue_RobustMr-Egger','BritishMR_Beta_RobustMr-Egger','BritishMR_SE_RobustMr-Egger',
            'BritishMR_Pvalue_SimpleMedian','BritishMR_Beta_SimpleMedian','BritishMR_SE_SimpleMedian',
            'BritishMR_Pvalue_WeightedMedian', 'BritishMR_Beta_WeightedMedian','BritishMR_SE_WeightedMedian',
            'BritishMR-IVWDelta_Model','BritishMR-IVWDelta_Pvalue','BritishMR-IVWDelta_Beta','BritishMR-IVWDelta_SE', 
            'BritishMR-IVWDelta_CILower','BritishMR-IVWDelta_CIUpper', 'BritishMR-IVWDelta_Alpha','BritishMR-IVWDelta_Robust', 'BritishMR-IVWDelta_Correlation',
            'BritishMR-IVWDelta_Penalized','BritishMR-IVWDelta_RSE','BritishMR-IVWDelta_Heter.Stat']

presso=['mr_presso_Pvalue','mr_presso_GlobalTest_Pvalue','mr_presso_Beta','mr_presso_Sd','mr_presso_T-stat','mr_presso_RSSobs']

heterogenity=['TwoSampleMR_heterogeneity_Q_Ivw','TwoSampleMR_heterogeneity_Q_df_Ivw','TwoSampleMR_heterogeneity_Q_pval_Ivw',      
              'TwoSampleMR_heterogeneity_Q_Ivw(Fixedeffects)', 'TwoSampleMR_heterogeneity_Q_df_Ivw(Fixedeffects)','TwoSampleMR_heterogeneity_Q_pval_Ivw(Fixedeffects)',
              'TwoSampleMR_heterogeneity_Q_Ivw(M-Randomeffects)','TwoSampleMR_heterogeneity_Q_df_Ivw(M-Randomeffects)','TwoSampleMR_heterogeneity_Q_pval_Ivw(M-Randomeffects)',
              'TwoSampleMR_heterogeneity_Q_Ivwradial','TwoSampleMR_heterogeneity_Q_df_Ivwradial', 'TwoSampleMR_heterogeneity_Q_pval_Ivwradial',
              'TwoSampleMR_heterogeneity_Q_Maximumlikelihood','TwoSampleMR_heterogeneity_Q_df_Maximumlikelihood','TwoSampleMR_heterogeneity_Q_pval_Maximumlikelihood',
              'TwoSampleMR_heterogeneity_Q_Mregger','TwoSampleMR_heterogeneity_Q_df_Mregger','TwoSampleMR_heterogeneity_Q_pval_Mregger',
              'TwoSampleMR_heterogeneity_Q_Mregger(Bootstrap)','TwoSampleMR_heterogeneity_Q_df_Mregger(Bootstrap)','TwoSampleMR_heterogeneity_Q_pval_Mregger(Bootstrap)',
              'TwoSampleMR_heterogeneity_Q_Penalisedweightedmedian', 'TwoSampleMR_heterogeneity_Q_df_Penalisedweightedmedian','TwoSampleMR_heterogeneity_Q_pval_Penalisedweightedmedian',
              'TwoSampleMR_heterogeneity_Q_Signconcordancetest', 'TwoSampleMR_heterogeneity_Q_df_Signconcordancetest','TwoSampleMR_heterogeneity_Q_pval_Signconcordancetest',
              'TwoSampleMR_heterogeneity_Q_Simplemedian','TwoSampleMR_heterogeneity_Q_df_Simplemedian','TwoSampleMR_heterogeneity_Q_pval_Simplemedian',
              'TwoSampleMR_heterogeneity_Q_Simplemode', 'TwoSampleMR_heterogeneity_Q_df_Simplemode','TwoSampleMR_heterogeneity_Q_pval_Simplemode',
              'TwoSampleMR_heterogeneity_Q_Simplemode(Nome)','TwoSampleMR_heterogeneity_Q_df_Simplemode(Nome)','TwoSampleMR_heterogeneity_Q_pval_Simplemode(Nome)',
              'TwoSampleMR_heterogeneity_Q_Unweightedregression','TwoSampleMR_heterogeneity_Q_df_Unweightedregression', 'TwoSampleMR_heterogeneity_Q_pval_Unweightedregression',
              'TwoSampleMR_heterogeneity_Q_Waldratio','TwoSampleMR_heterogeneity_Q_df_Waldratio','TwoSampleMR_heterogeneity_Q_pval_Waldratio',
              'TwoSampleMR_heterogeneity_Q_Weightedmedian','TwoSampleMR_heterogeneity_Q_df_Weightedmedian','TwoSampleMR_heterogeneity_Q_pval_Weightedmedian',
              'TwoSampleMR_heterogeneity_Q_Weightedmode','TwoSampleMR_heterogeneity_Q_df_Weightedmode','TwoSampleMR_heterogeneity_Q_pval_Weightedmode',
              'TwoSampleMR_heterogeneity_Q_Weightedmode(Nome)','TwoSampleMR_heterogeneity_Q_df_Weightedmode(Nome)', 'TwoSampleMR_heterogeneity_Q_pval_Weightedmode(Nome)']

hpleio=['TwoSampleMR_Hpleiotropy_egger_intercept', 'TwoSampleMR_Hpleiotropy_se', 'TwoSampleMR_Hpleiotropy_pval', 'TwoSampleMR_Direction_snp_r2.exposure', 
        'TwoSampleMR_Direction_snp_r2.outcome', 'TwoSampleMR_Direction_correct_causal_direction', 'TwoSampleMR_Direction_steiger_pval']
MR_Columns=TwoSampleMR+BritishMR+presso+heterogenity+hpleio



cis_table1a=pd.read_csv("/Users/jibinjohn/Downloads/UKBB-PPP_Decode_All_significannt_cis_exposure_AfterQC_LDclumping_Protein_level.csv")
cis_table1a[['HGNC','NCBI_Gene_ID']]=cis_table1a[['HGNC','NCBI_Gene_ID']].fillna(0).astype("int")
cis_table1a=cis_table1a.rename(columns={"exposure_UKBB-PPP":'UKBB_PPP_exposure',"exposure_Decode":'DECODE_exposure'})
cis_table1a=cis_table1a.fillna("NA")
cis_table1a['uniprot_ID']=cis_table1a['uniprot_ID'].str.split("-",expand=True)[0]
cis_ukb=cis_table1a[['uniprot_ID', 'Approved symbol', 'HGNC', 'Ensembl_gene_ID','NCBI_Gene_ID','id_UKBB-PPP', 'UKBB_PPP_exposure', 'oid_UKBB-PPP','platform_UKBB-PPP']].drop_duplicates()
cis_ukb=cis_ukb[cis_ukb['UKBB_PPP_exposure']!="NA"]
cis_decode=cis_table1a[['uniprot_ID', 'Approved symbol', 'HGNC', 'Ensembl_gene_ID','NCBI_Gene_ID','id_Decode', 'DECODE_exposure', 'platform_Decode','seqid_Decode']].drop_duplicates()
cis_decode=cis_decode[cis_decode['DECODE_exposure']!="NA"]
cis_ukb_decode=pd.merge(cis_ukb,cis_decode,on=['uniprot_ID', 'Approved symbol', 'HGNC', 'Ensembl_gene_ID','NCBI_Gene_ID'],how="outer")
cis_ukb_decode=pd.read_csv("/Users/jibinjohn/Downloads/UKBB-PPP_Decode_All_significannt_cis_exposure_AfterQC_LDclumping_Protein_level.csv")

trans_table1a=pd.read_csv("/Users/jibinjohn/Downloads/UKBB-PPP_Decode_All_significannt_trans_exposure_AfterQC_LDclumping_Protein_level.csv")
trans_table1a[['HGNC','NCBI_Gene_ID']]=trans_table1a[['HGNC','NCBI_Gene_ID']].fillna(0).astype("int")
trans_table1a=trans_table1a.rename(columns={"exposure_UKBB-PPP":'UKBB_PPP_exposure',"exposure_Decode":'DECODE_exposure'})
trans_table1a=trans_table1a.fillna("NA")
trans_table1a['uniprot_ID']=trans_table1a['uniprot_ID'].str.split("-",expand=True)[0]
trans_ukb=trans_table1a[['uniprot_ID', 'Approved symbol', 'HGNC', 'Ensembl_gene_ID','NCBI_Gene_ID','id_UKBB-PPP', 'UKBB_PPP_exposure', 'oid_UKBB-PPP','platform_UKBB-PPP']].drop_duplicates()
trans_ukb=trans_ukb[trans_ukb['UKBB_PPP_exposure']!="NA"]
trans_decode=trans_table1a[['uniprot_ID', 'Approved symbol', 'HGNC', 'Ensembl_gene_ID','NCBI_Gene_ID','id_Decode', 'DECODE_exposure', 'platform_Decode','seqid_Decode']].drop_duplicates()
trans_decode=trans_decode[trans_decode['DECODE_exposure']!="NA"]
trans_ukb_decode=pd.merge(trans_ukb,trans_decode,on=['uniprot_ID', 'Approved symbol', 'HGNC', 'Ensembl_gene_ID','NCBI_Gene_ID'],how="outer")
trans_ukb_decode=pd.read_csv("/Users/jibinjohn/Downloads/UKBB-PPP_Decode_All_significannt_trans_exposure_AfterQC_LDclumping_Protein_level.csv")

decode_platfor=pd.concat([cis_decode,trans_decode]).drop_duplicates()
decode_platfor_2=decode_platfor[['DECODE_exposure', 'platform_Decode','seqid_Decode']].drop_duplicates()
ukb_platfor=pd.concat([cis_ukb,trans_ukb]).drop_duplicates()
ukb_platfor_2=ukb_platfor[['UKBB_PPP_exposure','oid_UKBB-PPP','platform_UKBB-PPP']].drop_duplicates()



Results_df=pd.DataFrame()

for pqtltype in ['CisExposure','TransExposureNoMHC']:
    for disease in ["BIP_PGC3_noukb","Cognition","Depression_iPSYCH_2023","PGC3_SCZ_NoUKB"]:
        uk_scz=pd.read_csv(f"Biogen_{disease}_{pqtltype}_CompleteMR_AnalysisResults_withUniprotID.csv")
        de_scz=pd.read_csv(f"Decode_{disease}_{pqtltype}_CompleteMR_AnalysisResults_withUniprotID.csv")
        
        uk_scz_1=uk_scz[uk_common_columns]
        uk_scz_1=uk_scz_1.rename(columns={'BritishMR_SNPs':"UKB-PPP_BritishMR_SNPs",'BritishMR_nSNPs':"UKB-PPP_BritishMR_nSNPs",
                            'TwoSampleMR_nsnp_Ivw':"UKB-PPP_TwoSampleMR_nsnp_Ivw",'TwoSampleMR_nsnp_WaldRatio':"UKB-PPP_TwoSampleMR_nsnp_WaldRatio"})
        uk_scz_2=uk_scz[MR_Columns]
        uk_scz_2.columns=["UKB-PPP_"+x for x in uk_scz_2.columns ]
        uk_scz=pd.concat([uk_scz_1,uk_scz_2],axis=1)
        uk_scz=uk_scz.rename(columns={"exposure":"UKB-PPP_exposure"})
        
        decode_common_columns=["Approved symbol","DECODE_gene_name","UniProt_ID","HGNC","DECODE_exposure","Gene_Symbol","outcome","exposure", 'CIS/TRANS']
        de_scz_1=de_scz[decode_common_columns]
        de_scz_1=de_scz_1.rename(columns={'BritishMR_SNPs':"Decode_BritishMR_SNPs",'BritishMR_nSNPs':"Decode_BritishMR_nSNPs",
                            'TwoSampleMR_nsnp_Ivw':"Decode_TwoSampleMR_nsnp_Ivw",'TwoSampleMR_nsnp_WaldRatio':"Decode_TwoSampleMR_nsnp_WaldRatio"})
        de_scz_2=de_scz[MR_Columns]
        de_scz_2.columns=["Decode_"+x for x in de_scz_2.columns ]
        de_scz=pd.concat([de_scz_1,de_scz_2],axis=1)
        de_scz=de_scz.rename(columns={"exposure":"Decode_exposure"})
        
        merged_results=pd.merge(uk_scz,de_scz,
          on=['Approved symbol','UniProt_ID', 'HGNC', 'outcome', 'CIS/TRANS'],how="outer")
        Results_df=pd.concat([Results_df,merged_results])



def combinepvalue(row):
    p1=row['UKB-PPP_BritishMR-IVWDelta_Pvalue']
    p2=row['Decode_BritishMR-IVWDelta_Pvalue']
    result=combine_pvalues([p1, p2], method='fisher')[1]
    return result

Results_df1=Results_df[Results_df['UKB-PPP_BritishMR-IVWDelta_Pvalue'].notna() & Results_df['Decode_BritishMR-IVWDelta_Pvalue'].notna()]
Results_df2=Results_df[Results_df['UKB-PPP_BritishMR-IVWDelta_Pvalue'].isna() | Results_df['Decode_BritishMR-IVWDelta_Pvalue'].isna()]

Results_df1['metap_Pvalue'] = Results_df1.apply(combinepvalue, axis=1)
Results_df=pd.concat([Results_df1,Results_df2])

Results_df=pd.merge(Results_df,decode_platfor_2,on="DECODE_exposure",how="left")
Results_df=pd.merge(Results_df,ukb_platfor_2,on="UKBB_PPP_exposure",how="left")


dbm_cis_df=Results_df.copy()

conditions2 = [(dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'].isna() & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] >= 0)).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] >= 0) & dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'].isna()).astype(bool),
                (dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'].isna() & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] < 0)).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] < 0) & dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'].isna()).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] < 0) & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] < 0)).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] >= 0) & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] >= 0)).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] >= 0) & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] < 0)).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] < 0) & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] >= 0)).astype(bool)]
        
values2 = ["NA,+", "+,NA", "NA,-", "-,NA", "-,-", "+,+", "+,-", "-,+"]
dbm_cis_df['Estimate_Direction'] = np.select(conditions2, values2, default="Default_Value")

conditions = [((~dbm_cis_df['DECODE_exposure'].isna()) & (~dbm_cis_df['UKBB_PPP_exposure'].isna())),
            (dbm_cis_df['DECODE_exposure'].isna()) & (~dbm_cis_df['UKBB_PPP_exposure'].isna()),
            (~dbm_cis_df['DECODE_exposure'].isna()) & (dbm_cis_df['UKBB_PPP_exposure'].isna())]

values = ['BOTH', 'UKBB-PPP_Specific', 'DECODE_Specific']
dbm_cis_df['Cateegary'] = np.select(conditions, values)

dbm_cis_df["Lowest_Pvalue"]=dbm_cis_df[['Decode_BritishMR-IVWDelta_Pvalue', 'UKB-PPP_BritishMR-IVWDelta_Pvalue', 'metap_Pvalue']].apply(lambda row: row.min(), axis=1)
dbm_cis_df["Lowest_Pvalue_Source"]=dbm_cis_df[['Decode_BritishMR-IVWDelta_Pvalue', 'UKB-PPP_BritishMR-IVWDelta_Pvalue', 'metap_Pvalue']].idxmin(axis=1).str.replace("_Pvalue","")
dbm_cis_df["Largest_Pvalue"]=dbm_cis_df[['Decode_BritishMR-IVWDelta_Pvalue', 'UKB-PPP_BritishMR-IVWDelta_Pvalue', 'metap_Pvalue']].apply(lambda row: row.max(), axis=1)
dbm_cis_df["Largest_Pvalue_Source"]=dbm_cis_df[['Decode_BritishMR-IVWDelta_Pvalue', 'UKB-PPP_BritishMR-IVWDelta_Pvalue', 'metap_Pvalue']].idxmax(axis=1).str.replace("_Pvalue","")
dbm_cis_df["Lowest_Pvalue_Source"]=dbm_cis_df["Lowest_Pvalue_Source"].str.replace("_BritishMR-IVWDelta","")
dbm_cis_df["Largest_Pvalue_Source"]=dbm_cis_df["Largest_Pvalue_Source"].str.replace("_BritishMR-IVWDelta","")
dbm_cis_df=dbm_cis_df.rename(columns={"Gene_Symbol_x":"UKBB_PPP_Gene_Symbol","Gene_Symbol_x":"Decode_Gene_Symbol"})

dbm_cis_df.to_csv("TwoSampleMR_BritishMR_SCZ_BIP_MDD_CTP_Analaysis_AllMethod.csv",index=None)





Results_df2=dbm_cis_df.copy()
Results_df2['UKB-PPP_nsnp']=np.where(Results_df2['UKB-PPP_TwoSampleMR_nsnp_Ivw'].notna(),Results_df2['UKB-PPP_TwoSampleMR_nsnp_Ivw'],Results_df2['UKB-PPP_TwoSampleMR_nsnp_WaldRatio'])
Results_df2['Decode_nsnp']=np.where(Results_df2['Decode_TwoSampleMR_nsnp_Ivw'].notna(),Results_df2['Decode_TwoSampleMR_nsnp_Ivw'],Results_df2['Decode_TwoSampleMR_nsnp_WaldRatio'])

common_cols=["Approved symbol","UniProt_ID","HGNC","outcome",'CIS/TRANS','Cateegary','Estimate_Direction','Lowest_Pvalue', 'Lowest_Pvalue_Source', 'Largest_Pvalue','Largest_Pvalue_Source','metap_Pvalue']

ukb_snp_cols=['UKB-PPP_nsnp','UKB-PPP_BritishMR_SNPs']
decode_snp_cols=['Decode_nsnp','Decode_BritishMR_SNPs']

ukb_ivw_dlta_cols=['UKB-PPP_BritishMR-IVWDelta_Model', 'UKB-PPP_BritishMR-IVWDelta_Pvalue', 'UKB-PPP_BritishMR-IVWDelta_Beta', 'UKB-PPP_BritishMR-IVWDelta_SE', 
                    'UKB-PPP_BritishMR-IVWDelta_CILower', 'UKB-PPP_BritishMR-IVWDelta_CIUpper', 'UKB-PPP_BritishMR-IVWDelta_Alpha', 'UKB-PPP_BritishMR-IVWDelta_Robust', 
                    'UKB-PPP_BritishMR-IVWDelta_Correlation', 'UKB-PPP_BritishMR-IVWDelta_Penalized', 'UKB-PPP_BritishMR-IVWDelta_RSE', 'UKB-PPP_BritishMR-IVWDelta_Heter.Stat']

decode_ivw_dlta_cols=['Decode_BritishMR-IVWDelta_Model', 'Decode_BritishMR-IVWDelta_Pvalue', 'Decode_BritishMR-IVWDelta_Beta', 'Decode_BritishMR-IVWDelta_SE', 
                      'Decode_BritishMR-IVWDelta_CILower', 'Decode_BritishMR-IVWDelta_CIUpper', 'Decode_BritishMR-IVWDelta_Alpha', 
                      'Decode_BritishMR-IVWDelta_Robust', 'Decode_BritishMR-IVWDelta_Correlation', 'Decode_BritishMR-IVWDelta_Penalized', 'Decode_BritishMR-IVWDelta_RSE', 'Decode_BritishMR-IVWDelta_Heter.Stat']

two_smr_cols=['TwoSampleMR_Pvalue_Ivw','TwoSampleMR_Beta_Ivw', 'TwoSampleMR_SE_Ivw','TwoSampleMR_Pvalue_Ivw(FixedEffects)','TwoSampleMR_Beta_Ivw(FixedEffects)',
              'TwoSampleMR_SE_Ivw(FixedEffects)','TwoSampleMR_Pvalue_Ivw(M-Randomeffects)','TwoSampleMR_Beta_Ivw(M-Randomeffects)','TwoSampleMR_SE_Ivw(M-Randomeffects)',
              'TwoSampleMR_Pvalue_WaldRatio','TwoSampleMR_Beta_WaldRatio','TwoSampleMR_SE_WaldRatio']

het_cols=['TwoSampleMR_heterogeneity_Q_Ivw(Fixedeffects)','TwoSampleMR_heterogeneity_Q_df_Ivw(Fixedeffects)',
          'TwoSampleMR_heterogeneity_Q_pval_Ivw(Fixedeffects)','TwoSampleMR_heterogeneity_Q_Ivw(M-Randomeffects)',
          'TwoSampleMR_heterogeneity_Q_df_Ivw(M-Randomeffects)','TwoSampleMR_heterogeneity_Q_pval_Ivw(M-Randomeffects)',
          'TwoSampleMR_heterogeneity_Q_Waldratio','TwoSampleMR_heterogeneity_Q_df_Waldratio','TwoSampleMR_heterogeneity_Q_pval_Waldratio']

ukb_ivw_het_cols=["UKB-PPP_"+x for x in het_cols]
decode_het_cols=["Decode_"+x for x in het_cols]

ukb_ivw_hpleio_cols=["UKB-PPP_"+x for x in hpleio]
decode_hpleio_cols=["Decode_"+x for x in hpleio]

ukb_ivw_two_smr_cols=["UKB-PPP_"+x for x in two_smr_cols]
decode_two_smr_cols=["Decode_"+x for x in two_smr_cols]

last_columns=['platform_Decode','seqid_Decode','oid_UKBB-PPP','platform_UKBB-PPP','UKBB_PPP_exposure',"DECODE_exposure"]
selectd_cols=common_cols+ukb_ivw_dlta_cols+ukb_snp_cols+decode_ivw_dlta_cols+decode_snp_cols+ukb_ivw_het_cols+decode_het_cols+ukb_ivw_hpleio_cols+decode_hpleio_cols+ukb_ivw_two_smr_cols+decode_two_smr_cols+last_columns
Results_df2=Results_df2[selectd_cols]

ukb_rename_columns={'UKB-PPP_BritishMR-IVWDelta_Model':'UKB-PPP_Model',      'UKB-PPP_BritishMR-IVWDelta_Pvalue':'UKB-PPP_Pvalue',           'UKB-PPP_BritishMR-IVWDelta_Beta':'UKB-PPP_Estimate', 'UKB-PPP_BritishMR-IVWDelta_SE':'UKB-PPP_StdError', 
                    'UKB-PPP_BritishMR-IVWDelta_CILower':'UKB-PPP_CILower',  'UKB-PPP_BritishMR-IVWDelta_CIUpper':'UKB-PPP_CIUpper',         'UKB-PPP_BritishMR-IVWDelta_Alpha':'UKB-PPP_Alpha', 
                    'UKB-PPP_BritishMR-IVWDelta_Robust':'UKB-PPP_Robust',    'UKB-PPP_BritishMR-IVWDelta_Correlation':'UKB-PPP_Correlation', 'UKB-PPP_BritishMR-IVWDelta_Penalized':'UKB-PPP_Penalized', 
                    'UKB-PPP_BritishMR-IVWDelta_RSE':'UKB-PPP_RSE',          'UKB-PPP_BritishMR-IVWDelta_Heter.Stat':'UKB-PPP_Heter.Stat'}

decod_rename_columns={'DECODE_BritishMR-IVWDelta_Model':'DECODE_Model',     'DECODE_BritishMR-IVWDelta_Pvalue':'DECODE_Pvalue',           'DECODE_BritishMR-IVWDelta_Beta':'DECODE_Estimate', 'DECODE_BritishMR-IVWDelta_SE':'DECODE_StdError', 
                      'DECODE_BritishMR-IVWDelta_CILower':'DECODE_CILower', 'DECODE_BritishMR-IVWDelta_CIUpper':'DECODE_CIUpper',         'DECODE_BritishMR-IVWDelta_Alpha':'DECODE_Alpha', 
                      'DECODE_BritishMR-IVWDelta_Robust':'DECODE_Robust',   'DECODE_BritishMR-IVWDelta_Correlation':'DECODE_Correlation', 'DECODE_BritishMR-IVWDelta_Penalized':'DECODE_Penalized', 
                      'DECODE_BritishMR-IVWDelta_RSE':'DECODE_RSE',         'DECODE_BritishMR-IVWDelta_Heter.Stat':'DECODE_Heter.Stat'}


Results_df2.to_csv("TwoSampleMR_BritishMR_SCZ_BIP_MDD_CTP_Analaysis_SelectedMethods_SuppTable2.csv",index=None)










#########################-----------------------------------------------------------------######################################################################





















UKB_PPP_selectd_cols=['UKB-PPP_'+x for x in selectd_cols]
Decode_selectd_cols=['Decode_'+x for x in selectd_cols]

Results_df2=Results_df[common_cols+snp_cols+UKB_PPP_selectd_cols+Decode_selectd_cols]

ivw_dlta_cols_rennamee={'BritishMR-IVWDelta_Model':'DECODE_Model','BritishMR-IVWDelta_Pvalue':'DECODE_Pvalue','BritishMR-IVWDelta_Beta':'DECODE_Estimate','BritishMR-IVWDelta_SE':'DECODE_StdError', 
                        'BritishMR-IVWDelta_CILower':'DECODE_CILower','BritishMR-IVWDelta_CIUpper':'DECODE_CIUpper','BritishMR-IVWDelta_Penalized':'DECODE_Penalized',

                        'BritishMR-IVWDelta_Alpha':'DECODE_Alpha','BritishMR-IVWDelta_Robust':'DECODE_Robust', 'BritishMR-IVWDelta_Correlation':'DECODE_Correlation',
                        ,'BritishMR-IVWDelta_RSE':'DECODE_RSE','BritishMR-IVWDelta_Heter.Stat':'DECODE_Heter.Stat'}



Results_df2.to_csv("TwoSampleMR_BritishMR_SCZ_BIP_MDD_CTP_Analaysis_SelectedMethods.csv",index=None)













Results_df2=Results_df2.rename(columns={"outcome":"Outcome"})
Results_df2[['UKBB_PPP_exposure','DECODE_exposure','Outcome','CIS/TRANS','HGNC','UniProt_ID']].fillna("NA",inplace=True)
Results_df2[['UKBB_PPP_exposure','DECODE_exposure']]=Results_df2[['UKBB_PPP_exposure','DECODE_exposure']].fillna("NA")



s_table1s_table=pd.read_excel("/Users/jibinjohn/Downloads/Supplementary_Table2_MR.xlsx",skiprows=2)
s_table[['UKBB_PPP_exposure',"DECODE_exposure"]]=s_table[['UKBB_PPP_exposure',"DECODE_exposure"]].fillna("NA")
s_table[['UKBB_PPP_exposure','DECODE_exposure','Outcome','CIS/TRANS','HGNC','UniProt_ID']].fillna("NA",inplace=True)




Results_df2=pd.merge(Results_df2,decode_platfor_2,on="DECODE_exposure",how="left")
Results_df2=pd.merge(Results_df2,ukb_platfor_2,on="UKBB_PPP_exposure",how="left")
Results_df2.to_csv("TwoSampleMR_BritishMR_SCZ_BIP_MDD_CTP_Analaysis_SelectedMethods_SupTable2.csv",index=None)

dbm_cis_df=Results_df2.copy()
dbm_cis_df=dbm_cis_df.rename(columns={"sum":"metap_Pvalue"})

conditions2 = [(dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'].isna() & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] >= 0)).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] >= 0) & dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'].isna()).astype(bool),
                (dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'].isna() & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] < 0)).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] < 0) & dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'].isna()).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] < 0) & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] < 0)).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] >= 0) & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] >= 0)).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] >= 0) & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] < 0)).astype(bool),
                ((dbm_cis_df['Decode_BritishMR-IVWDelta_Beta'] < 0) & (dbm_cis_df['UKB-PPP_BritishMR-IVWDelta_Beta'] >= 0)).astype(bool)]
        
values2 = ["NA,+", "+,NA", "NA,-", "-,NA", "-,-", "+,+", "+,-", "-,+"]
dbm_cis_df['Estimate_Direction'] = np.select(conditions2, values2, default="Default_Value")

conditions = [((~dbm_cis_df['DECODE_exposure'].isna()) & (~dbm_cis_df['UKBB_PPP_exposure'].isna())),
            (dbm_cis_df['DECODE_exposure'].isna()) & (~dbm_cis_df['UKBB_PPP_exposure'].isna()),
            (~dbm_cis_df['DECODE_exposure'].isna()) & (dbm_cis_df['UKBB_PPP_exposure'].isna())]

values = ['BOTH', 'UKBB-PPP_Specific', 'DECODE_Specific']
dbm_cis_df['Cateegary'] = np.select(conditions, values)

dbm_cis_df["Lowest_Pvalue"]=dbm_cis_df[['Decode_BritishMR-IVWDelta_Pvalue', 'UKB-PPP_BritishMR-IVWDelta_Pvalue', 'metap_Pvalue']].apply(lambda row: row.min(), axis=1)
dbm_cis_df["Lowest_Pvalue_Source"]=dbm_cis_df[['Decode_BritishMR-IVWDelta_Pvalue', 'UKB-PPP_BritishMR-IVWDelta_Pvalue', 'metap_Pvalue']].idxmin(axis=1).str.replace("_Pvalue","")
dbm_cis_df["Largest_Pvalue"]=dbm_cis_df[['Decode_BritishMR-IVWDelta_Pvalue', 'UKB-PPP_BritishMR-IVWDelta_Pvalue', 'metap_Pvalue']].apply(lambda row: row.max(), axis=1)
dbm_cis_df["Largest_Pvalue_Source"]=dbm_cis_df[['Decode_BritishMR-IVWDelta_Pvalue', 'UKB-PPP_BritishMR-IVWDelta_Pvalue', 'metap_Pvalue']].idxmax(axis=1).str.replace("_Pvalue","")

dbm_cis_df.to_csv("TwoSampleMR_BritishMR_SCZ_BIP_MDD_CTP_Analaysis_SelectedMethods_SupTable2.csv",index=None)









['UKB-PPP_BritishMR-IVWDelta_Beta']
['Decode_BritishMR-IVWDelta_Beta']

s_table1=pd.concat([s_table1a,s_table1b]).drop_duplicates()
#s_table1.to_csv("/Users/jibinjohn/Downloads/UKBB-PPP_Decode_All_significannt_CIS_trans_exposure_AfterQC_LDclumping_Protein_level.csv")
#s_table1=s_table1.rename(columns={"exposure_UKBB-PPP":'UKBB_PPP_exposure',"exposure_Decode":'DECODE_exposure'})


s_table1_1=pd.read_csv("/Users/jibinjohn/Downloads/UKBB-PPP_Decode_All_significannt_CIS_trans_exposure_AfterQC_LDclumping_Protein_level.csv")
s_table1_1=s_table1_1.rename(columns={"exposure_UKBB-PPP":'UKBB_PPP_exposure',"exposure_Decode":'DECODE_exposure'})
s_table1_1[['HGNC','NCBI_Gene_ID']]=s_table1_1[['HGNC','NCBI_Gene_ID']].fillna(0).astype("int")
s_table1_1=s_table1_1.fillna("NA")

merged=pd.merge(s_table1_1,Results_df2,on=['UKBB_PPP_exposure','DECODE_exposure'],how="outer")
merged.to_csv("TwoSampleMR_BritishMR_SCZ_BIP_MDD_CTP_Analaysis_SelectedMethods_SupTable2.csv",index=None)
