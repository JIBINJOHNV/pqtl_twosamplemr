import pandas as pd
import numpy as np
import os,glob


basedir=os.getcwd()+"/"
pqtltype="Biogen"
gwasnames=['PGC_ADHD2022_iPSYCH_deCODE', 'BIP_PGC3_noukb', 'ASD_PGC', 'PGC3_SCZ', 'PGC_MDD_Depression']
cis_trans=['TransExposureNoMHC', 'TransExposure', 'CisExposure', 'TransExposureNoMHCUnique']

os.system("mkdir CombinedResultsfromAllBatches")
for gwasname in gwasnames:
    for cistran in cis_trans:
        file_prefix=f"{pqtltype}_{gwasname}_{cistran}"
        print(file_prefix)
        unique_files=glob.glob(f"Part0/{file_prefix}_*")
        unique_files=[x.split("/")[1] for x in unique_files ]
        unique_files=[x for x in unique_files if x.endswith("csv")]
        print(unique_files)
        os.system(f"mkdir -p CombinedResultsfromAllBatches/{gwasname}/{file_prefix}")
        file_cout=len(glob.glob(f"Part*/{file_prefix}_Harmonised_Exposure_Outcome.csv"))
        if file_cout==50:
            for file in unique_files:
                Files=glob.glob(f"Part*/{file}")
                master_df=pd.DataFrame()
                for part_file in Files:
                    part_file_df=pd.read_csv(part_file)
                    master_df=pd.concat([master_df,part_file_df])
                master_df.to_csv(f"CombinedResultsfromAllBatches/{gwasname}/{file_prefix}/{file}",index=None)
        else:
            print(f"Run not yet completed only {file_prefix}Harmonised_Exposure_Outcome.csv are {file_cout}")

               


import pandas as pd
import os,glob


basedir=os.getcwd()+"/"
pqtltype="Biogen"
gwasnames=['BIP_PGC3_noukb']
cis_trans=['TransExposureNoMHC', 'TransExposure', 'CisExposure', 'TransExposureNoMHCUnique']


for gwasname in gwasnames:
    for cistran in cis_trans:
        file_prefix=f"{pqtltype}_{gwasname}_{cistran}"
        #TwoSampleMR_Analysis_Multiple_MR_Test
        TwoSampleMultipleMRTest_file=glob.glob(f"CombinedResultsfromAllBatches/{gwasname}/{file_prefix}/*TwoSampleMR_Analysis_Multiple_MR_Test.csv")[0]
        TwoSampleMultipleMRTest_df=pd.read_csv(TwoSampleMultipleMRTest_file)
        TwoSampleMultipleMRTest_df=TwoSampleMultipleMRTest_df.drop(["id.exposure", "id.outcome"],axis=1).drop_duplicates()
        TwoSampleMultipleMRTest_df=TwoSampleMultipleMRTest_df.sort_values(["outcome", "exposure","method","pval"]).drop_duplicates(subset=["outcome", "exposure", "method"],keep='first')
        TwoSampleMultipleMRTest_df["method"]=TwoSampleMultipleMRTest_df["method"].str.replace("Inverse variance weighted","IVW").str.replace("multiplicative random effects","M-RandomEffects").str.title().str.replace(" ","")
        
        MR_Analysis=TwoSampleMultipleMRTest_df.pivot(index=["outcome" ,"exposure"],values=["pval","nsnp","b","se","lo_ci","up_ci","or","or_lci95","or_uci95"],columns="method").reset_index()
        MR_Analysis.columns=['_'.join(col) for col in MR_Analysis.columns.values]
        MR_Analysis.rename(columns={'outcome_':'outcome', 'exposure_':'exposure'},inplace=True)
        MR_Analysis.columns=['outcome', 'exposure']+["TwoSampleMR_"+x for x in MR_Analysis.columns if x not in ['outcome', 'exposure']]
        MR_Analysis_pvaluiecolumns=MR_Analysis[[x for x in MR_Analysis.columns if all(not x.startswith(y) for y in ["TwoSampleMR_lo_ci", "TwoSampleMR_up_ci", "TwoSampleMR_or", "TwoSampleMR_or_lci95", "TwoSampleMR_or_uci95","TwoSampleMR_se_","TwoSampleMR_b_"])]]
        nsnp_columns=[x for x in MR_Analysis_pvaluiecolumns.columns if x.startswith("TwoSampleMR_nsnp_")]
        MR_Analysis_pvaluiecolumns["TwoSampleMR_nsnps"]=MR_Analysis_pvaluiecolumns[nsnp_columns].max(axis=1)
        MR_Analysis_pvaluiecolumns=MR_Analysis_pvaluiecolumns.drop(nsnp_columns,axis=1)
        MR_Analysis_pvaluiecolumns.columns=['outcome', 'exposure']+["TwoSampleMR_"+x for x in MR_Analysis_pvaluiecolumns.columns if x not in ['outcome', 'exposure']]
        
        ####Hetrogenity
        het_file=glob.glob(f"CombinedResultsfromAllBatches/{gwasname}/{file_prefix}/*TwoSampleMR_Analysis_HeterogenityTest.csv")[0]
        het=pd.read_csv(het_file)
        het=het.drop(["id.exposure", "id.outcome"],axis=1)
        het=het.sort_values(["outcome", "exposure","method","Q_pval"]).drop_duplicates(subset=["outcome", "exposure", "method"],keep='first')
        het = het.drop_duplicates(subset=["outcome", "exposure", "method"],keep='first')
        het["method"]=TwoSampleMultipleMRTest_df["method"].str.replace("Inverse variance weighted","IVW").str.replace("multiplicative random effects","M-RandomEffects").str.title().str.replace(" ","")
        het=het[~het['exposure'].isna()]
        het=het.drop_duplicates(subset=["outcome","exposure","method"])
        
        het_2=het.pivot(index=["outcome" ,"exposure"],values=["Q","Q_df","Q_pval"],columns="method").reset_index()
        het_2.columns=['_'.join(col) for col in het_2.columns.values]
        het_2.columns=['outcome', 'exposure']+["heterogeneity_"+x.replace(" ","") for x in het_2.columns if x not in ['id.exposure_', 'id.outcome_', 'outcome_', 'exposure_']]
        het_2.columns=['outcome', 'exposure']+["TwoSampleMR_"+x for x in het_2.columns if x not in ['outcome', 'exposure']]
        het_2_pvaluiecolumns=het_2[['outcome', 'exposure']+[x for x in het_2.columns if x.startswith("TwoSampleMR_heterogeneity_Q_pval_")]]
        het_2_pvaluiecolumns.columns=['outcome', 'exposure']+["TwoSampleMR_"+x for x in het_2_pvaluiecolumns.columns if x not in ['outcome', 'exposure']]
        
        ###pleio
        hpleio_file=glob.glob(f"CombinedResultsfromAllBatches/{gwasname}/{file_prefix}/*TwoSampleMR_Analysis_Hpleiotropy_Test.csv")[0]
        pleio=pd.read_csv(hpleio_file)
        pleio=pleio[~pleio['exposure'].isna()]
        pleio=pleio.drop(["id.exposure", "id.outcome"],axis=1)
        pleio=pleio.sort_values(["outcome", "exposure","pval"])
        pleio = pleio.drop_duplicates(subset=["outcome", "exposure"],keep='first')
        pleio.columns=["outcome" ,"exposure"]+["TwoSampleMR_Hpleiotropy_"+x for x in pleio.columns if x not in ["outcome" ,"exposure"]]
        pleio_pvaluiecolumns=pleio[['outcome', 'exposure','TwoSampleMR_Hpleiotropy_pval']]
        
        ##Direction
        dir_file=glob.glob(f"CombinedResultsfromAllBatches/{gwasname}/{file_prefix}/*TwoSampleMR_Analysis_Directionality_Test.csv")[0]
        dirction=pd.read_csv(dir_file)
        dirction=dirction[~dirction['exposure'].isna()]
        dirction=dirction.drop(["id.exposure", "id.outcome"],axis=1)
        dirction=dirction.sort_values(["outcome", "exposure","steiger_pval"])
        dirction = dirction.drop_duplicates(subset=["outcome", "exposure"],keep='first')
        dirction.columns=["exposure" ,"outcome"]+["TwoSampleMR_Direction_"+x for x in dirction.columns if x not in ["outcome" ,"exposure"]]
        dirction_pvaluiecolumns=dirction[['outcome', 'exposure','TwoSampleMR_Direction_correct_causal_direction']]
        
        ##Merge All TwosampleMR Results
        mr_het=pd.merge(MR_Analysis,het_2,on=['outcome', 'exposure'],how="outer")
        mr_het_ple=pd.merge(mr_het,pleio,on=['outcome', 'exposure'],how="outer")
        mr_het_ple_dir=pd.merge(mr_het_ple,dirction,on=[ 'outcome', 'exposure'],how="outer")
        
        ##Merge All TwosampleMR Results with relevant columns (P values)
        mr_pleiohet_pvaluiecolumns=pd.merge(pleio_pvaluiecolumns,het_2_pvaluiecolumns,on=['outcome', 'exposure'],how="outer")
        mr_pleiohetdirecton_pvaluiecolumns=pd.merge(mr_pleiohet_pvaluiecolumns,dirction_pvaluiecolumns,on=['outcome', 'exposure'],how="outer")
        
        #### Biogen pipeline
        mr_file=glob.glob(f"CombinedResultsfromAllBatches/{gwasname}/{file_prefix}/*MendelianPipelineTest.csv")[0]
        mrpipelie=pd.read_csv(mr_file)
        mrpipelie=mrpipelie[~mrpipelie['exposure'].isna()]
        mrpipelie=mrpipelie.drop(["id.exposure", "id.outcome"],axis=1)
        mrpipelie=mrpipelie.sort_values(["outcome", "exposure",'pval'])
        inv=mrpipelie[mrpipelie['snp'].str.contains(",")]
        inv=inv.sort_values(by=["exposure","outcome",'pval'])
        inv = inv.drop_duplicates(subset=["outcome", "exposure"],keep='first')
        
        inv2=inv[['exposure', 'outcome']]
        inv2["type"]="INV"
        merged=pd.merge(mrpipelie,inv2,on=['exposure', 'outcome'],how="outer")
        wald=merged[merged["type"]!="INV"]
        mr_pipeline_result=pd.concat([wald,inv]).drop_duplicates()
        mr_pipeline_result.drop("type",axis=1,inplace=True)
        mr_pipeline_result.columns=["exposure" ,"outcome"]+["Biogen_"+x for x in mr_pipeline_result.columns if x not in ["outcome" ,"exposure"]]
        mr_pipeline_result_pvaluiecolumns=mr_pipeline_result[['exposure', 'outcome','Biogen_method', 'Biogen_nsnp','Biogen_pval']]
        
        ################################-------------------------- Mendelian randomisation all test
        MRAlltest_files=glob.glob(f"CombinedResultsfromAllBatches/{gwasname}/{file_prefix}/*MendelianRandomization_AllTest.csv")[0]
        MRAlltest_df=pd.read_csv(MRAlltest_files)
        MRAlltest_df=MRAlltest_df[~MRAlltest_df["Method"].str.contains("intercept")]
        MRAlltest_df.columns=['Method', 'Estimate', 'StdError', '95%CI_Lower','95%CI_Upper','P-value', 'outcome','exposure', 'SNPs']
        MRAlltest_df["Method"]=MRAlltest_df["Method"].str.title().str.replace(" ","")
        MRAlltest_Result_df=MRAlltest_df.pivot(index=["outcome" ,"exposure",'SNPs'],values=['Estimate', 'StdError', '95%CI_Lower','95%CI_Upper','P-value'],columns="Method").reset_index()
        MRAlltest_Result_df.columns=['_'.join(col) for col in MRAlltest_Result_df.columns.values]
        MRAlltest_Result_df.rename(columns={'outcome_':'outcome', 'exposure_':'exposure',"SNPs_":"SNPs"},inplace=True)
        MRAlltest_Result_df["nSNPs"]=MRAlltest_Result_df['SNPs'].str.count(',')+1
        MRAlltest_Result_df.columns=["outcome","exposure"]+["BritishMR_"+x for x in MRAlltest_Result_df.columns if x not in ["outcome" ,"exposure"]]
        MRAlltest_Result_df_pvaluiecolumns=MRAlltest_Result_df[['exposure', 'outcome', 'BritishMR_nSNPs']+[x for x in MRAlltest_Result_df.columns if "BritishMR_P-value_" in x ]]
        
        ################################-------------------------- Mendelian randomisation IVW test
        MRIVWtest_files=glob.glob(f"CombinedResultsfromAllBatches/{gwasname}/{file_prefix}/*MendelianRandomization_IVW_Delta_Test.csv")[0]
        MRIVWtest_df=pd.read_csv(MRIVWtest_files)
        MRIVWtest_df=MRIVWtest_df[['Exposure','Outcome','Model','Penalized', 'Estimate', 'CILower', 'Alpha', 'SNPs', 'Heter.Stat','Robust', 'Correlation', 'StdError', 'CIUpper', 'Pvalue', 'RSE']]
        MRIVWtest_df.columns=['Exposure','Outcome']+["BritishMR-IVWDelta_"+x for x in MRIVWtest_df.columns if x not in ['Outcome','Exposure']]
        MRIVWtest_df.rename(columns={'Outcome':'outcome','Exposure':'exposure'},inplace=True)
        MRIVWtest_df_pvaluiecolumns=MRIVWtest_df[['outcome', 'exposure','BritishMR-IVWDelta_SNPs','BritishMR-IVWDelta_Heter.Stat','BritishMR-IVWDelta_Pvalue']]
        
        #########################---------------------------------Presso----------------------------------------------------- 
        presso_files=glob.glob(f"CombinedResultsfromAllBatches/{gwasname}/{file_prefix}/*MendelianRandomization_Presso_Test.csv")[0]
        MRPRESSO_df=pd.read_csv(presso_files)
        MRPRESSO_raw_df=MRPRESSO_df[MRPRESSO_df["MR Analysis"]=="Raw"].drop(["Exposure", "MR Analysis"],axis=1)
        MRPRESSO_raw_df.columns=["mr_presso_"+x for x in MRPRESSO_raw_df.columns if x not in ['exposure', 'outcome'] ]+['exposure', 'outcome']
        MRPRESSO_raw_df_pvaluiecolumns=MRPRESSO_raw_df[['outcome', 'exposure','mr_presso_P-value']]
        
        #########Merge With All columns 
        twosamplemr_MRPRESSO=pd.merge(mr_het_ple_dir,MRPRESSO_raw_df,on=["outcome" ,"exposure"],how="outer")
        biogen_british_ivwdelta=pd.merge(mr_pipeline_result,MRIVWtest_df,on=["outcome" ,"exposure"],how="outer")
        britishall_biogen_britishivwdelta=pd.merge(MRAlltest_Result_df,biogen_british_ivwdelta,on=["outcome", "exposure"],how="outer")
        britishall_biogen_britishivwdelta_twosamplemr_MRPRESSO=pd.merge(twosamplemr_MRPRESSO,britishall_biogen_britishivwdelta,on=["outcome", "exposure"],how="outer")
        colorder=list(britishall_biogen_britishivwdelta_twosamplemr_MRPRESSO.columns)
        britishall_biogen_britishivwdelta_twosamplemr_MRPRESSO["Gene_Symbol"]=britishall_biogen_britishivwdelta_twosamplemr_MRPRESSO["exposure"].str.split("_",expand=True)[0]
        britishall_biogen_britishivwdelta_twosamplemr_MRPRESSO=britishall_biogen_britishivwdelta_twosamplemr_MRPRESSO[['Gene_Symbol']+colorder]
        britishall_biogen_britishivwdelta_twosamplemr_MRPRESSO.to_csv(f"{file_prefix}_CompleteMR_AnalysisResults.csv",index=None)
        
        ##Merge only selected important columns
        mr_pleiohet_pvaluiecolumns=pd.merge(pleio_pvaluiecolumns,het_2_pvaluiecolumns,on=['outcome', 'exposure'],how="outer")
        mr_pleiohetdirecton_pvaluiecolumns=pd.merge(mr_pleiohet_pvaluiecolumns,dirction_pvaluiecolumns,on=['outcome', 'exposure'],how="outer")
        biogen_british_ivwDelta=pd.merge(mr_pipeline_result_pvaluiecolumns,MRIVWtest_df_pvaluiecolumns,on=["exposure","outcome"],how="outer")
        biogen_british_ivwDelta_allbritish=pd.merge(biogen_british_ivwDelta,MRAlltest_Result_df_pvaluiecolumns,on=["exposure","outcome"],how="outer")
        biogen_british_ivwDelta_allbritish_presso=pd.merge(biogen_british_ivwDelta_allbritish,MRPRESSO_raw_df_pvaluiecolumns,on=["exposure","outcome"],how="outer")
        allmr_pvalues=pd.merge(biogen_british_ivwDelta_allbritish_presso,MR_Analysis_pvaluiecolumns,on=['outcome', 'exposure'],how="outer")
        no_pvaluecolumns=['Biogen_nsnp','BritishMR-IVWDelta_SNPs','BritishMR_nSNPs','BritishMR-IVWDelta_Heter.Stat','Biogen_method','TwoSampleMR_TwoSampleMR_nsnps']
        allmr_pvalues_pvaluesonly=allmr_pvalues.drop(no_pvaluecolumns,axis=1)
        allmr_pvalues_nopvalues=allmr_pvalues[no_pvaluecolumns]
        allmr_pvalues=pd.concat([allmr_pvalues_pvaluesonly,allmr_pvalues_nopvalues],axis=1)
        allmr_pvalues_dirhetpleio=pd.merge(allmr_pvalues,mr_pleiohetdirecton_pvaluiecolumns,on=['outcome', 'exposure'],how="outer")
        allmr_pvalues_dirhetpleio.columns=[ x.replace("TwoSampleMR_TwoSampleMR_","TwoSampleMR_") for x in allmr_pvalues_dirhetpleio.columns ]
        colorder=list(allmr_pvalues_dirhetpleio.columns)
        allmr_pvalues_dirhetpleio["Gene_Symbol"]=allmr_pvalues_dirhetpleio["exposure"].str.split("_",expand=True)[0]
        allmr_pvalues_dirhetpleio=allmr_pvalues_dirhetpleio[['Gene_Symbol']+colorder]
        allmr_pvalues_dirhetpleio.to_csv(f"{file_prefix}_CompleteMR_AnalysisResults_WithselectedColumns.csv",index=None)
