import pandas as pd 


gwases = ["ASD_PGC", "BIP_PGC3_noukb", "Depression_iPSYCH_2023", "PGC3_SCZ", "PGC_ADHD2022_iPSYCH_deCODE", "PGC_AN2"]
pqtl_types = ["Biogen", "Decode"]
exposure_types=["TransExposureNoMHC","CisExposure"]

for gwas in gwases:
    final_result_df=pd.DataFrame(columns=['Gene_Symbol', 'SNP'])
    for exposure_type in exposure_types:
        result_df = pd.DataFrame(columns=['Gene_Symbol', 'SNP'])
        for pqtl_type in pqtl_types:
            file = f"{pqtl_type}_{gwas}_{exposure_type}_CompleteMR_AnalysisResults.xlsx"
            df = pd.read_excel(file, sheet_name="Singlevariant")
            df = df[['Gene_Symbol', 'SNP', 'outcome', 'TwosampleMR-wald_Pvalue', 'TwosampleMR-wald_Beta']]
            df = df.sort_values(by=['TwosampleMR-wald_Pvalue', 'SNP', 'Gene_Symbol']).drop_duplicates(
                subset=["Gene_Symbol", "SNP"], keep="first")
            df = df.sort_values(by="SNP")
            count_df = df.groupby("SNP")["outcome"].count().reset_index().rename(columns={"outcome": "SNP_Count"})
            df = pd.merge(df, count_df, on="SNP")
            df.columns = ['Gene_Symbol', "SNP"] + [x + f"_{gwas}_{pqtl_type}" for x in df.columns if x not in ['Gene_Symbol', "SNP"]]
            har_df = pd.read_excel(file, sheet_name="harmonised_data")
            har_df=har_df[[ 'SNP', 'gene.exposure', 'eaf.exposure','eaf.outcome']].sort_values(by=["SNP","gene.exposure","eaf.exposure"])
            har_df=har_df.drop_duplicates(subset=["SNP","gene.exposure"],keep="first")
            har_df=har_df.rename(columns={'gene.exposure': 'Gene_Symbol'})
            har_df.columns = ["SNP",'Gene_Symbol'] + [x + f"_{gwas}_{pqtl_type}" for x in har_df.columns if  x not in ['Gene_Symbol', "SNP"]]
            df=pd.merge(df,har_df,on=['Gene_Symbol', 'SNP'], how="left")
            result_df = pd.merge(result_df, df, on=['Gene_Symbol', 'SNP'], how="outer")
        first_part_columns = ['Gene_Symbol', 'SNP'] + [x for x in result_df.columns if "SNP_Count" in x]
        second_part_columns = [x for x in result_df.columns if "wald_Pvalue" in x]
        third_part_columns = [x for x in result_df.columns if x not in first_part_columns + second_part_columns]
        result_df = result_df[first_part_columns + second_part_columns + third_part_columns]
        count_columns=[x for x in result_df.columns if "SNP_Count_" in x]
        pvalue_columns=[x for x in result_df.columns if 'TwosampleMR-wald_Pvalu' in x]
        result_df=result_df[(result_df[pvalue_columns[0]]<=0.001) | (result_df[pvalue_columns[1]]<=0.001) ]
        result_df=result_df[(result_df[count_columns[0]]>1) | (result_df[count_columns[1]]>1) ]
        result_df.to_csv(f"{gwas}_{exposure_type}_NonUniqVariants_Pvalue_0.001.csv",index=None)
        fial_column_names=[x.replace('TwosampleMR-wald_',"") for x in result_df.columns]
        result_df.columns=fial_column_names
        result_df.columns=['Gene_Symbol',"SNP"] + [f"{exposure_type}_"+x for x in result_df.columns if  x not in ['Gene_Symbol', "SNP"]]
        final_result_df=pd.merge(final_result_df,result_df,on=['Gene_Symbol', 'SNP'], how="outer")
        first_columns = ['Gene_Symbol', 'SNP'] + [x for x in final_result_df.columns if "SNP_Count" in x]
        second_columns = [x for x in final_result_df.columns if "wald_Pvalue" in x]
        third_columns = [x for x in final_result_df.columns if "_Pvalue_" in x]
        forth_columns = [x for x in final_result_df.columns if "_Beta_" in x]
        fifth_columns = [x for x in final_result_df.columns if "eaf.exposure" in x]
        sixth_columns = [x for x in final_result_df.columns if "eaf.outcome" in x]
        seventh_columns = [x for x in final_result_df.columns if x not in first_columns + second_columns +third_columns+forth_columns+fifth_columns+sixth_columns]
        final_result_df = final_result_df[first_columns + second_columns +third_columns+forth_columns+fifth_columns+sixth_columns+seventh_columns]
        final_result_df.to_csv(f"{gwas}_NonUniqVariants_Pvalue_0.001.csv",index=None)
