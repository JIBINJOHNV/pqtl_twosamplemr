import pandas as pd
import os,glob

Prefix="Biogen_trans_exposure_noMHC_"  # "Biogen_CisExposure_" "Biogen_trans_exposure_NoMHC_Unique_" "Biogen_trans_exposure_"

cild=glob.glob("*TwoSampleMR_Analysis_Multiple_MR_Test.csv")
MR_Analysis=pd.DataFrame()

for file in cild:
     tempd_df=pd.read_csv(file)
     MR_Analysis=pd.concat([MR_Analysis,tempd_df])

MR_Analysis=MR_Analysis.drop(["id.exposure", "id.outcome"],axis=1)
MR_Analysis=MR_Analysis.drop_duplicates()
MR_Analysis=MR_Analysis.sort_values(["outcome", "exposure","method","pval"])
MR_Analysis = MR_Analysis.drop_duplicates(subset=["outcome", "exposure", "method"],keep='first')

MR_Analysis=MR_Analysis.pivot(index=["outcome" ,"exposure"],values=["pval","nsnp","b","se","lo_ci","up_ci","or","or_lci95","or_uci95"],columns="method").reset_index()
MR_Analysis.columns=['_'.join(col) for col in MR_Analysis.columns.values]
MR_Analysis.columns=['outcome', 'exposure']+[x.replace(" ","") for x in MR_Analysis.columns if x not in ['id.exposure_', 'id.outcome_', 'outcome_', 'exposure_']]


####Hetrogenity
het_file=glob.glob("*TwoSampleMR_Analysis_HeterogenityTest.csv")
het=pd.DataFrame()

for file in het_file:
     tempd_df=pd.read_csv(file)
     het=pd.concat([het,tempd_df])

het=het.drop(["id.exposure", "id.outcome"],axis=1)
het=het.drop_duplicates()
het=het.sort_values(["outcome", "exposure","method","Q_pval"])
het = het.drop_duplicates(subset=["outcome", "exposure", "method"],keep='first')

het=het[~het['exposure'].isna()]
het_2=het.pivot(index=["outcome" ,"exposure"],values=["Q","Q_df","Q_pval"],columns="method").reset_index()
het_2.columns=['_'.join(col) for col in het_2.columns.values]
het_2.columns=['outcome', 'exposure']+["heterogeneity_"+x.replace(" ","") for x in het_2.columns if x not in ['id.exposure_', 'id.outcome_', 'outcome_', 'exposure_']]


###pleio
hpleio_file=glob.glob("*TwoSampleMR_Analysis_Hpleiotropy_Test.csv")
pleio=pd.DataFrame()

for file in hpleio_file:
     tempd_df=pd.read_csv(file)
     pleio=pd.concat([pleio,tempd_df])

pleio=pleio[~pleio['exposure'].isna()]
pleio=pleio.drop(["id.exposure", "id.outcome"],axis=1)

pleio=pleio.sort_values(["outcome", "exposure","pval"])
pleio = pleio.drop_duplicates(subset=["outcome", "exposure"],keep='first')
pleio.columns=["outcome" ,"exposure"]+["Hpleiotropy_"+x for x in pleio.columns if x not in ["outcome" ,"exposure"]]


##Direction
dir_file=glob.glob("*TwoSampleMR_Analysis_Directionality_Test.csv")
dirction=pd.DataFrame()

for file in dir_file:
     tempd_df=pd.read_csv(file)
     dirction=pd.concat([dirction,tempd_df])

dirction=dirction[~dirction['exposure'].isna()]
dirction=dirction.drop(["id.exposure", "id.outcome"],axis=1)
dirction=dirction.sort_values(["outcome", "exposure","steiger_pval"])
dirction = dirction.drop_duplicates(subset=["outcome", "exposure"],keep='first')
dirction.columns=["exposure" ,"outcome"]+["Direction_"+x for x in dirction.columns if x not in ["outcome" ,"exposure"]]

mr_het=pd.merge(MR_Analysis,het_2,on=['outcome', 'exposure'],how="outer")
mr_het_ple=pd.merge(mr_het,pleio,on=['outcome', 'exposure'],how="outer")
mr_het_ple_dir=pd.merge(mr_het_ple,dirction,on=[ 'outcome', 'exposure'],how="outer")


#### mIRPIPELINE
mr_file=glob.glob("*MendelianPipelineTest.csv")
mrpipelie=pd.DataFrame()

for file in mr_file:
     tempd_df=pd.read_csv(file)
     mrpipelie=pd.concat([mrpipelie,tempd_df])

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

mr_pipeline_result.columns=["exposure" ,"outcome"]+["MR_Pipeline_"+x for x in mr_pipeline_result.columns if x not in ["outcome" ,"exposure"]]
mr_pipeline_result.drop("MR_Pipeline_type",axis=1,inplace=True)


####################################------------------------------------Mendelian randomisation-----------------------------------------------------###########################################

################################-------------------------- Mendelian randomisation all test
MRAlltest_files=glob.glob("*MendelianRandomization_AllTest.csv")
MRAlltest_df=pd.DataFrame()

for file in MRAlltest_files:
     tempd_df=pd.read_csv(file)
     MRAlltest_df=pd.concat([MRAlltest_df,tempd_df])

MRAlltest_df.columns=['Method', 'Estimate', 'Std Error', '95% CI_Lower','95% CI_Upper','P-value', 'outcome','exposure', 'SNPs']
MRAlltest_Result_df=MRAlltest_df[["exposure","outcome","SNPs"]].drop_duplicates()

for method in [x for x in MRAlltest_df["Method"].unique() if "intercept" not in x]:
     mrtemdf=MRAlltest_df[MRAlltest_df["Method"]==method].drop("Method",axis=1)
     mrtemdf.columns=[method.replace(" ","")+"_"+x for x in mrtemdf.columns if x not in ['outcome', 'exposure', 'SNPs']]+['outcome', 'exposure', 'SNPs']
     MRAlltest_Result_df=pd.merge(MRAlltest_Result_df,mrtemdf,on=['outcome', 'exposure', 'SNPs'],how="outer")

MRAlltest_Result_df["nSNPs"]=MRAlltest_Result_df['SNPs'].str.count(',')+1
MRAlltest_Result_df.rename(columns={"exposure":"Exposure","outcome":"Outcome"},inplace=True)
                           
################################-------------------------- Mendelian randomisation IVW test
MRIVWtest_files=glob.glob("*MendelianRandomization_IVW_Delta_Test.csv")
MRIVWtest_df=pd.DataFrame()

for file in MRIVWtest_files:
     tempd_df=pd.read_csv(file)
     MRIVWtest_df=pd.concat([MRIVWtest_df,tempd_df])

MRIVWtest_df=MRIVWtest_df[['Outcome','Exposure','Estimate','StdError','Pvalue','CILower','CIUpper','SNPs']]
MRIVWtest_df.columns=['Outcome','Exposure']+["MRIVWtest_"+x for x in MRIVWtest_df.columns if x not in ['Outcome','Exposure']]




#########################---------------------------------Presso----------------------------------------------------- 
presso_files=glob.glob("*MendelianRandomization_Presso_Test.csv")
MRPRESSO_df=pd.DataFrame()

for file in presso_files:
     tempd_df=pd.read_csv(file)
     MRPRESSO_df=pd.concat([MRPRESSO_df,tempd_df])

MRPRESSO_raw_df=MRPRESSO_df[MRPRESSO_df["MR Analysis"]=="Raw"].drop(["Exposure", "MR Analysis"],axis=1)
MRPRESSO_raw_df.columns=[ x for x in MRPRESSO_raw_df.columns if x not in ['exposure', 'outcome'] ]+['exposure', 'outcome']
MRPRESSO_raw_df.rename(columns={"exposure":"Exposure","outcome":"Outcome"},inplace=True)
MRPRESSO_raw_df.columns=["mr_presso_"+x for x in MRPRESSO_raw_df.columns if x not in ['Outcome','Exposure']]+['Exposure','Outcome']



#########Merge 

mr_pipeline_mrResults=pd.merge(mr_het_ple_dir,mr_pipeline_result,on=["outcome" ,"exposure"],how="outer")
mr_pipeline_mrResults.rename(columns={"exposure":"Exposure","outcome":"Outcome"},inplace=True)


MR_All_df=pd.merge(MRAlltest_Result_df,MRIVWtest_df,on=["Outcome", "Exposure"],how="outer")
MR_All_df_MRPRESSO_raw_df=pd.merge(MR_All_df,MRPRESSO_raw_df,on=['Exposure','Outcome'],how="outer")

Final_mrResults=pd.merge(mr_pipeline_mrResults,MR_All_df_MRPRESSO_raw_df,on=["Outcome", "Exposure"],how="outer")
Final_mrResults["Gene_Symbol"]=Final_mrResults["Exposure"].str.split("_",expand=True)[0]
Final_mrResults.to_csv(f"{Prefix}CompleteMR_AnalysisResults.csv",index=None)



# List of strings to search for
patterns_to_check = ["MRIVWtest_", "MR_Pipeline", "Hpleiotropy", "heterogeneity", "Direction"]
Final_mrResults2=Final_mrResults[["Outcome","Exposure","Gene_Symbol"]+[ x for x in Final_mrResults.columns if [y for y in patterns_to_check if y in x ]  ]]
Final_mrResults2.columns=["Outcome","Exposure","Gene_Symbol"]+[ x for x in Final_mrResults2.columns if  not [y for y in ["Outcome","Exposure","Gene_Symbol"] if y in x ]  ]
Final_mrResults2.to_csv(f"{Prefix}CompleteMR_AnalysisResults_ForMtaP.csv",index=None)




#################----------------------------------------Single variant---------------------------------------------------
##Metafixed
cild=glob.glob("*TwoSampleMR_Analysis_SingleVariantMetafixed_Test.csv")
meta=pd.DataFrame()

for file in cild:
     tempd_df=pd.read_csv(file)
     meta=pd.concat([meta,tempd_df])

meta=meta[meta["SNP"].str.contains("^rs")]
meta=meta.drop(["id.exposure", "id.outcome"],axis=1)
meta=meta.drop_duplicates()
meta=meta.sort_values(["outcome", "exposure","SNP","p"])
meta = meta.drop_duplicates(subset=["outcome", "exposure","SNP"],keep='first')
meta=meta.rename(columns={"b":"metafixed_b","se":"metafixed_se","p":"metafixed_p"})

##Wald test
cild=glob.glob("*TwoSampleMR_Analysis_SingleVariantWald_Test.csv")
wald=pd.DataFrame()

for file in cild:
     tempd_df=pd.read_csv(file)
     wald=pd.concat([wald,tempd_df])

wald=wald[wald["SNP"].str.contains("^rs")]
wald=wald.drop(["id.exposure", "id.outcome"],axis=1)
wald=wald.drop_duplicates()
wald=wald.sort_values(["outcome", "exposure","SNP","p"])
wald = wald.drop_duplicates(subset=["outcome", "exposure","SNP"],keep='first')
wald=wald.rename(columns={"b":"wald_b","se":"wald_se","p":"wald_p"})

##Mr pipeline
mr_file=glob.glob("*MendelianPipelineTest.csv")
mrpipelie=pd.DataFrame()

for file in mr_file:
     tempd_df=pd.read_csv(file)
     mrpipelie=pd.concat([mrpipelie,tempd_df])

mrpipelie=mrpipelie[~mrpipelie['exposure'].isna()]
mrpipelie=mrpipelie.drop(["id.exposure", "id.outcome"],axis=1)
mrpipelie2=mrpipelie[mrpipelie["method"]=="Wald ratio"]
mrpipelie2=mrpipelie2.sort_values(["exposure","outcome","method","nsnp","snp",'pval'])
mrpipelie2=mrpipelie2.drop_duplicates(subset=["exposure","outcome","method","nsnp","snp"],keep="first")

mrpipelie2.columns=["exposure" ,"outcome"]+["MR_Pipeline_"+x for x in mrpipelie2.columns if x not in ["outcome" ,"exposure"]]
meta_wald=pd.merge(meta,wald,on=['exposure', 'outcome', 'samplesize', 'SNP'])

singlevariant_mr=pd.merge(meta_wald,mrpipelie2,left_on=['exposure','outcome','SNP'],right_on=["exposure","outcome","MR_Pipeline_snp"],how="outer")
singlevariant_mr["Gene_Symbol"]=singlevariant_mr["exposure"].str.split(":",expand=True)[0]

singlevariant_mr.to_csv(f"{Prefix}SingleVariant_CompleteMR_AnalysisResults.csv",index=None)
