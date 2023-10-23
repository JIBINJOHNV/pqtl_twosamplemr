import pandas as pd
import numpy as np


protein_info_df=pd.read_csv("Complete_UKBB-PPP_Deecode_protein_withcomplete_Geneinformation_SeelecteedColumns_Annotation_Drug_cellstr_size_striped.csv")
decode_df=protein_info_df[['DECODE_exposure','UniProt_ID']].drop_duplicates()
decode_df=decode_df[~decode_df['DECODE_exposure'].isna()]
biogen_df=protein_info_df[['UKBB_PPP_exposure','UniProt_ID']].drop_duplicates()
biogen_df=biogen_df[~biogen_df['UKBB_PPP_exposure'].isna()]



PqtlSource=["Decode","Biogen"]
gwases=["Depression_iPSYCH_2023" ,"PGC3_SCZ_NoUKB" ,"Cognition" , "BIP_PGC3_noukb","PGC_ADHD2022_iPSYCH_deCODE","PGC_AN2","ASD_PGC"]
cis_trans_decode_biogen_df=pd.DataFrame()
pqtls=["CisExposure","TransExposureNoMHC"]

        #gwas=gwases[0]          
        #pqtl=pqtls[0]   
result_df=pd.DataFrame()

for gwas in gwases:
    for pqtl in pqtls:
        biomet_df=pd.read_csv(f"Biogen_{gwas}_{pqtl}_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta.csv")
        biomet_df=pd.merge(biomet_df,biogen_df,left_on="exposure",right_on="UKBB_PPP_exposure",how="left").drop("UKBB_PPP_exposure",axis=1)
        biomet_df.to_csv(f"Biogen_{gwas}_{pqtl}_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta.csv",index=None)
        
        decomet_df=pd.read_csv(f"Decode_{gwas}_{pqtl}_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta.csv")
        decomet_df=pd.merge(decomet_df,decode_df,left_on="exposure",right_on="DECODE_exposure",how="left").drop("DECODE_exposure",axis=1)
        decomet_df.to_csv(f"Decode_{gwas}_{pqtl}_British_IVDeltaMRPipeline_AnalysisResults_WithselectedColumns_ForMeta.csv",index=None)

