import pandas as pd 
import numpy as np

df=pd.read_csv("Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2_withGeneInfo.csv")
#df=df[["Approved symbol","DECODE_gene_name","UKBB_PPP_gene_name","UniProt_ID","HGNC"]].drop_duplicates()
ukb_df=df[["Approved symbol","UKBB_PPP_Protein","UniProt_ID","HGNC","UKBB_PPP_exposure"]].drop_duplicates() ## UKBB_PPP_Protein UKBB_PPP_gene_name  ,"UKBB_PPP_exposure"
ukb_df=ukb_df[ukb_df["UKBB_PPP_Protein"].notna()].drop_duplicates()
data = {"Approved symbol": ["LILRA3"],"UKBB_PPP_Protein": ["LILRA3"],
        "UniProt_ID": ["Q8N6C8"],"HGNC": [6604], "UKBB_PPP_exposure":["LILRA3:Q8N6C8:Cardiometabolic_II"]}
ukb_df=pd.concat([ukb_df,pd.DataFrame(data)]).reset_index(drop=True)

df=pd.read_csv("Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2_withGeneInfo.csv")
decode_df=df[["Approved symbol","DECODE_gene_name","UniProt_ID","HGNC", 'DECODE_exposure']].drop_duplicates()
decode_df=decode_df[decode_df["DECODE_gene_name"].notna()].drop_duplicates()




ivw_wald_columns=['Approved symbol', 'DECODE_gene_name', 'UniProt_ID', 'HGNC', 'Gene_Symbol', 'outcome', 'exposure', 'CIS/TRANS', 
             'TwoSampleMR_Pvalue_Ivw(FixedEffects)', 'TwoSampleMR_nsnp_Ivw(FixedEffects)', 'TwoSampleMR_Beta_Ivw(FixedEffects)', 'TwoSampleMR_SE_Ivw(FixedEffects)','TwoSampleMR_Pvalue_WaldRatio', 'TwoSampleMR_nsnp_WaldRatio', 'TwoSampleMR_Beta_WaldRatio', 'TwoSampleMR_SE_WaldRatio']

phenotypes=['BIP_PGC3_noukb','Cognition','Depression_iPSYCH_2023','PGC3_SCZ_NoUKB']


for exposure_source in ["Decode", "Biogen"]:
    for pqtltype in ["TransExposureNoMHC", "CisExposure"]:
        for phenotype in phenotypes:
            df2 = pd.read_excel(f'{exposure_source}_{phenotype}_{pqtltype}_CompleteMR_AnalysisResults.xlsx', sheet_name="Complete_MR_Results")
            df2 = df2.drop("Unnamed: 0", axis=1)
            file = f'{exposure_source}_{phenotype}_{pqtltype}_CompleteMR_AnalysisResults.xlsx'
            out_name = file.replace(".xlsx", "")
            
            if exposure_source == "Biogen":
                a_symbol = pd.merge(ukb_df, df2, right_on=["Gene_Symbol",'exposure'], left_on=["UKBB_PPP_Protein","UKBB_PPP_exposure"], how="right")
            elif exposure_source == "Decode":
                a_symbol = pd.merge(decode_df, df2, right_on=["Gene_Symbol",'exposure'], left_on=["DECODE_gene_name","DECODE_exposure"], how="right")

            missing = a_symbol[a_symbol["Approved symbol"].isna()]

            if pqtltype == "CisExposure":
                a_symbol['CIS/TRANS'] = "CIS"
            elif pqtltype == "TransExposureNoMHC":
                a_symbol['CIS/TRANS'] = "TRANS"

            if missing.shape[0] == 0:
                not_to_include = ["_lo_ci_", "_up_ci_", "_or_", "_95%CI_", "_Biogen_up_ci", "_Biogen_lo_", "_Biogen_or_"]
                a_symbol = a_symbol[[x for x in a_symbol.columns if all(y not in x for y in not_to_include)]]
                a_symbol.columns = [x.replace(f"{phenotype}_", "") for x in a_symbol.columns]
                a_symbol.to_csv(f"{out_name}_withUniprotID.csv", index=None)
                a_symbol2 = a_symbol[[x for x in a_symbol.columns if any(x.endswith(y) for y in ivw_wald_columns)]]
                a_symbol2 = a_symbol2[[x for x in a_symbol2.columns if "_snp_r2." not in x]]
                a_symbol2.columns = [x.replace(f"TwoSampleMR_", "") for x in a_symbol2.columns]

                fixed_columns = ['Approved symbol', 'UniProt_ID', 'HGNC', 'Gene_Symbol', 'outcome', 'exposure', 'CIS/TRANS']
                ivw_columns = [x for x in a_symbol2.columns if "FixedEffects" in x]
                wald_columns = [x for x in a_symbol2.columns if "WaldRatio" in x]
                wald_df = a_symbol2[fixed_columns + wald_columns]
                wald_df.columns = [x.replace("_WaldRatio", "") for x in wald_df.columns]
                wald_df = wald_df[wald_df['Pvalue'].notna()]
                wald_df["method"] = "Wald ratio"
                ivw_df = a_symbol2[fixed_columns + ivw_columns]
                ivw_df.columns = [x.replace("_Ivw(FixedEffects)", "") for x in ivw_df.columns]
                ivw_df = ivw_df[ivw_df['Pvalue'].notna()]
                ivw_df["method"] = "Ivw(FixedEffects)"
                wald_ivw_df = pd.concat([wald_df, ivw_df], ignore_index=True)

                if exposure_source == "Decode":
                    Prefix = "Decode"
                    rename_column2 = {"Pvalue": f"{Prefix}_Pvalue", "nsnp": f"{Prefix}_nsnp",
                                        f"Beta": f"{Prefix}_Beta", "SE": f"{Prefix}_SE",
                                        "method": f"{Prefix}_method", "outcome": "Outcome",
                                        "exposure": "DECODE_exposure"}
                elif exposure_source == "Biogen":
                    Prefix = "UKB-PPP"
                    rename_column2 = {"Pvalue": f"{Prefix}_Pvalue", "nsnp": f"{Prefix}_nsnp",
                                        f"Beta": f"{Prefix}_Beta", "SE": f"{Prefix}_SE",
                                        "method": f"{Prefix}_method", "outcome": "Outcome",
                                        "exposure": "UKBB_PPP_exposure"}

                wald_ivw_df = wald_ivw_df.rename(columns=rename_column2)
                wald_ivw_df.to_csv(f"{out_name}_withUniprotID_IVW_Wald.csv", index=None)
            
            if missing.shape[0] != 0:
                print(f"Missing rows present in {out_name}",missing.shape[0] )




import pandas as pd 
import numpy as np

phenotypes=['BIP_PGC3_noukb','Cognition','Depression_iPSYCH_2023','PGC3_SCZ_NoUKB']

Final_df=pd.DataFrame()

for pqtltype in ["TransExposureNoMHC","CisExposure"]:
    for pqtltype in ["TransExposureNoMHC","CisExposure"]:
        for phenotype in phenotypes:
            decode_df=pd.read_csv(f'Decode_{phenotype}_{pqtltype}_CompleteMR_AnalysisResults_withUniprotID_IVW_Wald.csv')
            ukb_df=pd.read_csv(f'Biogen_{phenotype}_{pqtltype}_CompleteMR_AnalysisResults_withUniprotID_IVW_Wald.csv')
            merged=pd.merge(ukb_df,decode_df,on=['Approved symbol', 'UniProt_ID', 'HGNC', 'Outcome'],how="outer")
            merged=merged.drop(['CIS/TRANS_x','CIS/TRANS_y'],axis=1)
            if pqtltype=="CisExposure":
                merged['CIS/TRANS']="CIS"
            elif  pqtltype=="TransExposureNoMHC":
                merged['CIS/TRANS']="TRANS"
            
            merged=merged.rename(columns={"Gene_Symbol_x":"Decode_Gene_Symbol","Gene_Symbol_y":"UKB-PPP_Gene_Symbol"})
            merged[["Approved symbol","UniProt_ID","HGNC","DECODE_exposure","UKBB_PPP_exposure"]].fillna("NA",inplace=True)
            df=pd.read_csv("Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2_withGeneInfo.csv")
            df[["Approved symbol","UniProt_ID","HGNC","DECODE_exposure","UKBB_PPP_exposure"]].fillna("NA",inplace=True)
            merged=pd.merge(df,merged,on=["Approved symbol","UniProt_ID","HGNC","DECODE_exposure","UKBB_PPP_exposure","CIS/TRANS","Outcome"],how="right")
            merged.to_csv(f'UKB-PPP_Decode_{phenotype}_{pqtltype}_CompleteMR_AnalysisResults_withUniprotID_IVW_Wald.csv',index=None)
            Final_df=pd.concat([Final_df,merged]).drop_duplicates()


Final_df.to_csv(f'UKB-PPP_Decode_All_FoourPhenotype_CIS_TRANS_CompleteMR_AnalysisResults_withUniprotID_IVW_Wald.csv',index=None)

selected_columns=['Approved symbol', 'UniProt_ID', 'HGNC', 'CIS/TRANS', 'Outcome', 'Cateegary', 'Estimate_Direction', 'metap_Pvalue', 'Lowest_Pvalue', 'Largest_Pvalue', 'Lowest_Pvalue_Source', 'Largest_Pvalue_Source', 'DECODE_Model', 'DECODE_Pvalue', 'DECODE_Estimate', 'DECODE_StdError', 'DECODE_CILower', 'DECODE_CIUpper', 'DECODE_Penalized', 'DECODE_Robust', 'DECODE_Alpha', 'DECODE_Heter.Stat', 'DECODE_Correlation', 'DECODE_RSE', 'DECODE_N_SNPs', 'DECODE_SNPs', 'DECODE_snp_r2.exposure', 'DECODE_snp_r2.outcome', 'DECODE_correct_causal_direction', 'DECODE_steiger_pval', 'DECODE_egger_intercept', 'DECODE_se', 'DECODE_pval', 'DECODE_Q IVW', 'DECODE_Q IVW_radial', 'DECODE_Q MR_Egger', 'DECODE_Q Maximum_likelihood', 'DECODE_Q Unweighted_regression', 'DECODE_Q_df IVW', 'DECODE_Q_df IVW_radial', 'DECODE_Q_df MR_Egger', 'DECODE_Q_df Maximum_likelihood', 'DECODE_Q_df Unweighted_regression', 'DECODE_Q_pval IVW', 'DECODE_Q_pval IVW_radial', 'DECODE_Q_pval MR_Egger', 'DECODE_Q_pval Maximum_likelihood', 'DECODE_Q_pval Unweighted_regression', 'UKBB_PPP_Model', 'UKBB_PPP_Pvalue', 'UKBB_PPP_Estimate', 'UKBB_PPP_StdError', 'UKBB_PPP_CILower', 'UKBB_PPP_CIUpper', 'UKBB_PPP_Penalized', 'UKBB_PPP_Robust', 'UKBB_PPP_Alpha', 'UKBB_PPP_Heter.Stat', 'UKBB_PPP_Correlation', 'UKBB_PPP_RSE', 'UKBB_PPP_N_SNPs', 'UKBB_PPP_SNPs', 'UKBB_PPP_snp_r2.exposure', 'UKBB_PPP_snp_r2.outcome', 'UKBB_PPP_correct_causal_direction', 'UKBB_PPP_steiger_pval', 'UKBB_PPP_egger_intercept', 'UKBB_PPP_se', 'UKBB_PPP_pval', 'UKBB_PPP_Q IVW', 'UKBB_PPP_Q IVW_radial', 'UKBB_PPP_Q MR_Egger', 'UKBB_PPP_Q Maximum_likelihood', 'UKBB_PPP_Q Unweighted_regression', 'UKBB_PPP_Q_df IVW', 'UKBB_PPP_Q_df IVW_radial', 'UKBB_PPP_Q_df MR_Egger', 'UKBB_PPP_Q_df Maximum_likelihood', 'UKBB_PPP_Q_df Unweighted_regression', 'UKBB_PPP_Q_pval IVW', 'UKBB_PPP_Q_pval IVW_radial', 'UKBB_PPP_Q_pval MR_Egger', 'UKBB_PPP_Q_pval Maximum_likelihood', 'UKBB_PPP_Q_pval Unweighted_regression', 'UKB-PPP_Pvalue', 'UKB-PPP_nsnp', 'UKB-PPP_Beta', 'UKB-PPP_SE', 'UKB-PPP_method', 'UKB-PPP_Gene_Symbol', 'Decode_Pvalue', 'Decode_nsnp', 'Decode_Beta', 'Decode_SE', 'Decode_method', 'DECODE_seqid', 'DECODE_platform', 'DECODE_target_full_name', 'DECODE_target_name', 'UKBB_PPP_oid', 'UKBB_PPP_Panel', 'UKBB_PPP_platform', 'UKBB_PPP_Protein', 'UKBB_PPP_target_full_name', 'OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID', 'DECODE_exposure', 'UKBB_PPP_exposure']

Final_slected_columns_df=Final_df[selected_columns]
Final_slected_columns_df.to_csv(f'UKB-PPP_Decode_All_FoourPhenotype_CIS_TRANS_CompleteMR_AnalysisResults_withUniprotID_IVW_Wald_SelectedColumns.csv',index=None)
