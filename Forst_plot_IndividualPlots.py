
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import forestplot as fp


#https://github.com/LSYS/forestplot

pvalue_cutof=0.01
decimal_places = 2

# file location :/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/Forestplots/New_ForestPlots

df=pd.read_csv("Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2_withGeneInfo.csv")
#df=df[["Approved symbol","DECODE_gene_name","UKBB_PPP_gene_name","UniProt_ID","HGNC"]].drop_duplicates()
ukb_df=df[["Approved symbol","UKBB_PPP_Protein","UKBB_PPP_exposure"]].drop_duplicates() ## UKBB_PPP_Protein UKBB_PPP_gene_name  ,"UKBB_PPP_exposure"
ukb_df=ukb_df[ukb_df["UKBB_PPP_Protein"].notna()].drop_duplicates()
data = {"Approved symbol": ["LILRA3"],"UKBB_PPP_Protein": ["LILRA3"],"UKBB_PPP_exposure":["LILRA3:Q8N6C8:Cardiometabolic_II"]}
ukb_df=pd.concat([ukb_df,pd.DataFrame(data)]).reset_index(drop=True)

df=pd.read_csv("Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2_withGeneInfo.csv")
decode_df=df[["Approved symbol","DECODE_gene_name",'DECODE_exposure']].drop_duplicates()
decode_df=decode_df[decode_df["DECODE_gene_name"].notna()].drop_duplicates()




def save_mpl_fig(savename):
    savepath = os.path.join('/Users/jibinjohn/Desktop/MR_Analysis_Results_with_SCZ_NoUKB_Cognition_sup_table/Forestplots/New_ForestPlots/', savename)
    plt.savefig(f'{savepath}.png', dpi='figure', bbox_inches='tight')

def to_scientific(value):
    return '{:.{}e}'.format(value, decimal_places)


def  combined_mr(pqtl_type,excell,exposure_name):
    mr_df=pd.read_excel(excell,sheet_name="biogen_british")
    if pqtl_type=="decode":
        mr_df=mr_df[mr_df['exposure']==exposure_name]
    selected_columns=['Gene_Symbol','exposure','outcome','BritishMR-IVWDelta_Pvalue','BritishMR-IVWDelta_Beta','BritishMR-IVWDelta_SE','BritishMR-IVWDelta_nSNPs']
    include_columns=[x for  x in mr_df.columns if [y for y in selected_columns if y in x] ]
    mr_df2=mr_df[include_columns]

    pvalue_column=[x for x in include_columns if x.endswith("_Pvalue")][0]
    beta_column=[x for x in include_columns if x.endswith("_Beta") ][0]
    se_column=[x for x in include_columns if x.endswith("_SE") ][0]
    nsnp_columns=[x for x in include_columns if x.endswith("BritishMR-IVWDelta_nSNPs") ][0]
    mr_df2.rename(columns={pvalue_column:"P",beta_column:"b",se_column:"se","Gene_Symbol":"id.exposure",'outcome':'id.outcome',nsnp_columns:"nSNP" },inplace=True)

    if pqtl_type=="biogen":
        mr_df2['SNP']="Biogen"
        mr_df2["group"]="combined"
        mr_df2 = pd.merge(ukb_df, mr_df2, right_on=["id.exposure",'exposure'], left_on=["UKBB_PPP_Protein","UKBB_PPP_exposure"], how="right")
        mr_df2["id.exposure"]=mr_df2["Approved symbol"]
    
    elif pqtl_type=="decode":
        mr_df2['SNP']="Decode"
        mr_df2["group"]="combined"
        mr_df2 = pd.merge(decode_df, mr_df2, right_on=["id.exposure",'exposure'], left_on=["DECODE_gene_name","DECODE_exposure"], how="right")
        mr_df2["id.exposure"]=mr_df2["Approved symbol"]
    
    mr_df2['up'] = mr_df2['b'] + 1.96 * mr_df2['se']
    mr_df2['lo'] = mr_df2['b'] - 1.96 * mr_df2['se']
    mr_df2=mr_df2.sort_values(by=["id.exposure","P"])
    mr_df2=mr_df2.drop_duplicates(subset=["id.exposure"], keep="first")
    mr_df2["MAF"]="_"
    return mr_df2

def  single_mr(pqtl_type,excell,exposure_name):
    mr_df=pd.read_excel(excell,sheet_name="Singlevariant")
    if pqtl_type=="decode":
        mr_df=mr_df[mr_df['exposure']==exposure_name]
    selected_columns=['Gene_Symbol','exposure','outcome','TwosampleMR-wald_Beta', 'TwosampleMR-wald_SE','TwosampleMR-wald_Pvalue',"SNP"]
    mr_df2=mr_df[selected_columns]
    mr_df2.columns=['id.exposure','exposure', 'id.outcome','b', 'se','P',"SNP"]
    mr_df2['up'] = mr_df2['b'] + 1.96 * mr_df2['se']
    mr_df2['lo'] = mr_df2['b'] - 1.96 * mr_df2['se']
        
    if pqtl_type=="biogen":
        mr_df2["group"]="UKB-PPP_SingleSNP"
        mr_df2["SNP"]="UKB-PPP_"+mr_df2["SNP"]
        mr_df2 = pd.merge(ukb_df, mr_df2, right_on=['id.exposure','exposure'], left_on=["UKBB_PPP_Protein","UKBB_PPP_exposure"], how="right")
        mr_df2['id.exposure']=mr_df2["Approved symbol"]
                
    elif pqtl_type=="decode":
        mr_df2["group"]="Decode_SingleSNP"
        mr_df2["SNP"]="Decode_"+mr_df2["SNP"]
        mr_df2 = pd.merge(decode_df, mr_df2, right_on=['id.exposure','exposure'], left_on=["DECODE_gene_name","DECODE_exposure"], how="right")
        mr_df2['id.exposure']=mr_df2["Approved symbol"]
    
    mr_df2=mr_df2.sort_values(by=['id.exposure',"SNP"])
    mr_df2['nSNP']="_"
    
    if pqtl_type=="biogen":
        har_df=pd.read_excel(excell,sheet_name="harmonised_data")
        har_df["SNP"]="UKB-PPP_"+har_df["SNP"]
        har_df=har_df[["SNP","eaf.exposure","gene.exposure"]].drop_duplicates()
        har_df.rename(columns={"eaf.exposure":"MAF","gene.exposure":"id.exposure"},inplace=True)
        har_df["MAF"]=np.where(har_df["MAF"]<0.5,har_df["MAF"],1-har_df["MAF"])
        mr_df2=pd.merge(mr_df2,har_df,on=["SNP","id.exposure"])
    elif pqtl_type=="decode":
        har_df=pd.read_excel(excell,sheet_name="harmonised_data")
        har_df=har_df[har_df['exposure']==exposure_name]
        har_df=har_df[["SNP","eaf.exposure","gene.exposure"]].drop_duplicates()
        har_df.rename(columns={"eaf.exposure":"MAF","gene.exposure":"id.exposure"},inplace=True)
        har_df["MAF"]=np.where(har_df["MAF"]<0.5,har_df["MAF"],1-har_df["MAF"])
        har_df["SNP"]="Decode_"+har_df["SNP"]
        mr_df2=pd.merge(mr_df2,har_df,on=["SNP","id.exposure"])
    mr_df2["MAF"]=mr_df2["MAF"].round(4)
    return mr_df2



for pqtltype in ["CisExposure","TransExposureNoMHC"]:
    metap="Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2_withGeneInfo.csv"
    
    for Disease in ["Depression_iPSYCH_2023","BIP_PGC3_noukb","PGC3_SCZ_NoUKB","Cognition"]:
        Disease_dict={'Depression_iPSYCH_2023': 'Depression_iPSYCH_2023', 'BIP_PGC3_noukb': 'BIP_PGC3_noukb', 'PGC3_SCZ_NoUKB': 'PGC3_SCZ_NoUKB', 'Cognition': 'Cognition_Meta'}



Disease_dict={'Depression_iPSYCH_2023': 'Depression_iPSYCH_2023', 'BIP_PGC3_noukb': 'BIP_PGC3_noukb', 'PGC3_SCZ_NoUKB': 'PGC3_SCZ_NoUKB', 'Cognition': 'Cognition_Meta'}
metap="Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2_withGeneInfo.csv"        
pqtltype="TransExposureNoMHC"
Disease="PGC3_SCZ_NoUKB"

exposure_name="C15ORF48_6406_3"
exposure="C15ORF48"

os.system(f"mkdir -p {pqtltype}_BiogenDecode_Pqtl_{Disease}")
biogen_excell=f"Biogen_{Disease}_{pqtltype}_CompleteMR_AnalysisResults.xlsx"
decode_exsell=f"Decode_{Disease}_{pqtltype}_CompleteMR_AnalysisResults.xlsx"   
biogen_mr_df=combined_mr("biogen",biogen_excell,exposure_name)
decode_mr_df=combined_mr("decode",decode_exsell,exposure_name)
biogen_single_df=single_mr("biogen",biogen_excell,exposure_name)
decode_single_df=single_mr("decode",decode_exsell,exposure_name)
mr_df=pd.concat([biogen_mr_df,decode_mr_df])
mr_single=pd.concat([biogen_single_df,decode_single_df])
metap_df=pd.read_csv(metap)

if pqtltype=="TransExposureNoMHC" :
    metap_df=metap_df[metap_df['CIS/TRANS']=="TRANS"]
    metap_df=metap_df[metap_df['Outcome']==Disease_dict[Disease]]
    metap_df=metap_df[metap_df[f'Lowest_Pvalue']<pvalue_cutof]
    metap_df=metap_df[~metap_df[f'Lowest_Pvalue'].isna()]
    Exposures=list(set(metap_df[metap_df[f'Lowest_Pvalue']<pvalue_cutof]["Approved symbol"].unique()))
    Exposures.sort()

if pqtltype=="CisExposure" :
    metap_df=metap_df[metap_df['CIS/TRANS']=="CIS"]
    metap_df=metap_df[metap_df['Outcome']==Disease_dict[Disease]]
    metap_df=metap_df[metap_df[f'Lowest_Pvalue']<pvalue_cutof]
    metap_df=metap_df[~metap_df[f'Lowest_Pvalue'].isna()]
    Exposures=list(set(metap_df[metap_df[f'Lowest_Pvalue']<pvalue_cutof]["Approved symbol"].unique()))
    Exposures.sort()

ngenes=len(Exposures)
print(f" now running is {Disease} and {pqtltype} total  number of genes are {ngenes}")
single_selectd=mr_single[mr_single['id.exposure']==exposure].sort_values(by="SNP")
mr2_selectd=mr_df[mr_df['id.exposure']==exposure]
selected_df=pd.concat([mr2_selectd,single_selectd])
selected_df['up']=selected_df['up'].round(3)
selected_df['lo']=selected_df['lo'].round(3)
selected_df['b']=selected_df['b'].round(3)
selected_df['se']=selected_df['se'].round(3)
selected_df['P']=selected_df['P'].apply(to_scientific)

pqtls=list(mr2_selectd["SNP"].unique())
if len(pqtls)>1:
    out_name1="Both_"+"-".join(pqtls)
elif len(pqtls)==1:
    out_name1="Only_"+"-".join(pqtls)

lowest_pvalue_columns=[x for x in metap_df.columns if "Lowest_Pvalue" in x][0]
lowest_pvalue=metap_df[metap_df["Approved symbol"]==exposure][lowest_pvalue_columns].apply(to_scientific).unique()[0]
outname_final=out_name1+"_"+str(lowest_pvalue)

if selected_df.shape[0]<=29:
    fp.forestplot(selected_df, estimate="b",moerror="se",ll="lo", hl="up", varlabel="SNP",decimal_precision=3,ci_report=False,color_alt_rows=True,flush=False,table=True,groupvar="group",rightannote=['P', 'b','up', 'lo',"nSNP","MAF"],right_annoteheaders=["Pvalue", "Beta", "Up_CI","Low_CI","nSNP","MAF"],sort=True,**{"marker": "D", "markersize": 35,  "xlinestyle": (0, (10, 5)), "xlinecolor": "#808080", "xtick_size": 12}   )
    plt.savefig(f"{pqtltype}_BiogenDecode_Pqtl_{Disease}/{pqtltype}_BiogenDecode_Pqtl_{Disease}_{exposure}_plot_{outname_final}.png", bbox_inches="tight")

if selected_df.shape[0]>29:
    try:
        biogen=selected_df[selected_df["SNP"].str.contains("Biogen")]
        fp.forestplot(biogen, estimate="b",moerror="se",ll="lo", hl="up",varlabel="SNP",decimal_precision=3,
                ci_report=False,color_alt_rows=True,flush=False,table=True,groupvar="group",
                rightannote=['P', 'b','up', 'lo',"nSNP","MAF"],
                right_annoteheaders=["Pvalue", "Beta", "Up_CI","Low_CI","nSNP","MAF"],
                sort=True,**{"marker": "D", "markersize": 35,  "xlinestyle": (0, (10, 5)), "xlinecolor": "#808080", "xtick_size": 12}   )
        plt.savefig(f"{pqtltype}_BiogenDecode_Pqtl_{Disease}/{pqtltype}_BiogenDecode_Pqtl_{Disease}_{exposure}_BiogenPqtl_plot_{outname_final}.png", bbox_inches="tight")
    except:
        pass
    
    try:
        decode=selected_df[selected_df["SNP"].str.contains("Decode")]
        fp.forestplot(decode, estimate="b",moerror="se",ll="lo", hl="up",varlabel="SNP",decimal_precision=3,
                ci_report=False,color_alt_rows=True,flush=False,table=True,groupvar="group",
                rightannote=['P', 'b','up', 'lo',"nSNP","MAF"],
                right_annoteheaders=["Pvalue", "Beta", "Up_CI","Low_CI","nSNP","MAF"],
                sort=True,**{"marker": "D", "markersize": 35,  "xlinestyle": (0, (10, 5)), "xlinecolor": "#808080", "xtick_size": 12}   )
        plt.savefig(f"{pqtltype}_BiogenDecode_Pqtl_{Disease}/{pqtltype}_BiogenDecode_Pqtl_{Disease}_{exposure}_DecodePqtl_plot_{outname_final}.png", bbox_inches="tight")
    except:
        pass
