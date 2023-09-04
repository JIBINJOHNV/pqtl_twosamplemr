import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import forestplot as fp


#https://github.com/LSYS/forestplot

pvalue_cutof=0.01
decimal_places = 2

def save_mpl_fig(savename):
    savepath = os.path.join('/Users/jibinjohn/Downloads/BiogenUKB_Decode_Pqtl_MR/Exlsx/TransExposureNoMHC', savename)
    plt.savefig(f'{savepath}.png', dpi='figure', bbox_inches='tight')

def to_scientific(value):
    return '{:.{}e}'.format(value, decimal_places)


def  combined_mr(pqtl_type,excell):
    mr_df=pd.read_excel(excell,sheet_name="biogen_british")
    selected_columns=['Gene_Symbol','outcome','BritishMR-IVWDelta_Pvalue','BritishMR-IVWDelta_Beta','BritishMR-IVWDelta_SE','BritishMR-IVWDelta_nSNPs']
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
    
    elif pqtl_type=="decode":
        mr_df2['SNP']="Decode"
        mr_df2["group"]="combined"
    
    mr_df2['up'] = mr_df2['b'] + 1.96 * mr_df2['se']
    mr_df2['lo'] = mr_df2['b'] - 1.96 * mr_df2['se']
    mr_df2=mr_df2.sort_values(by=["id.exposure","P"])
    mr_df2=mr_df2.drop_duplicates(subset=["id.exposure"], keep="first")
    mr_df2["EAF"]="_"
    return mr_df2

def  single_mr(pqtl_type,excell):
    mr_df=pd.read_excel(excell,sheet_name="Singlevariant")
    selected_columns=['Gene_Symbol','outcome','TwosampleMR-wald_Beta', 'TwosampleMR-wald_SE','TwosampleMR-wald_Pvalue',"SNP"]
    mr_df2=mr_df[selected_columns]
    mr_df2.columns=['id.exposure', 'id.outcome','b', 'se','P',"SNP"]
    mr_df2['up'] = mr_df2['b'] + 1.96 * mr_df2['se']
    mr_df2['lo'] = mr_df2['b'] - 1.96 * mr_df2['se']
    
    if pqtl_type=="biogen":
        mr_df2["group"]="Biogen_SingleSNP"
        mr_df2["SNP"]="Biogen_"+mr_df2["SNP"]
    elif pqtl_type=="decode":
        mr_df2["group"]="Decode_SingleSNP"
        mr_df2["SNP"]="Decode_"+mr_df2["SNP"]
    mr_df2=mr_df2.sort_values(by=['id.exposure',"SNP"])
    mr_df2['nSNP']="_"
    
    har_df=pd.read_excel(excell,sheet_name="harmonised_data")
    har_df.rename(columns={"gene.exposure":"id.exposure"},inplace=True)
    har_df=har_df[["SNP","eaf.exposure","id.exposure"]].drop_duplicates()
    har_df.rename(columns={"eaf.exposure":"EAF"},inplace=True)
    if pqtl_type=="biogen":
        har_df["SNP"]="Biogen_"+har_df["SNP"]
    elif pqtl_type=="decode":
        har_df["SNP"]="Decode_"+har_df["SNP"]
    mr_df2=pd.merge(mr_df2,har_df,on=["SNP","id.exposure"])
    return mr_df2



for pqtltype in ["CisExposure","TransExposureNoMHC"]:
    metap="CIS_TransExposureNoMHC_Exposure_CompleteMR.xlsx"
    
    for Disease in ["PGC_AN2","Depression_iPSYCH_2023","BIP_PGC3_noukb","PGC3_SCZ" ,"ASD_PGC","PGC_ADHD2022_iPSYCH_deCODE"]:
        os.system(f"mkdir -p {pqtltype}_BiogenDecode_Pqtl_{Disease}")
        biogen_excell=f"Biogen_{Disease}_{pqtltype}_CompleteMR_AnalysisResults.xlsx"
        decode_exsell=f"Decode_{Disease}_{pqtltype}_CompleteMR_AnalysisResults.xlsx"   
        
        biogen_mr_df=combined_mr("biogen",biogen_excell)
        decode_mr_df=combined_mr("decode",decode_exsell)
        
        biogen_single_df=single_mr("biogen",biogen_excell)
        decode_single_df=single_mr("decode",decode_exsell)
        
        mr_df=pd.concat([biogen_mr_df,decode_mr_df])
        mr_single=pd.concat([biogen_single_df,decode_single_df])
        
        metap_df=pd.read_excel(metap)
        if Disease=="PGC_ADHD2022_iPSYCH_deCODE":
            metap_df.columns=[x.replace("ADHD2022_iPSYCH_deCODE",Disease) for x in metap_df.columns]
        
        gwas_columns=["Gene_Symbol"]+[ x for x in metap_df.columns if Disease in x]
        metap_df=metap_df[gwas_columns]
        
        if pqtltype=="TransExposureNoMHC" :
            gwas_columns=["Gene_Symbol"]+[ x for x in metap_df.columns if "TransExposureNoMHC_" in x]
            metap_df=metap_df[gwas_columns]
            metap_df=metap_df[metap_df[f'TransExposureNoMHC_{Disease}_LwestPvalue']<pvalue_cutof]
            metap_df=metap_df[~metap_df[f'TransExposureNoMHC_{Disease}_LwestPvalue'].isna()]
            Exposures=list(set(metap_df[metap_df[f'TransExposureNoMHC_{Disease}_LwestPvalue']<pvalue_cutof]["Gene_Symbol"].unique()))
        
        if pqtltype=="CisExposure" :
            gwas_columns=["Gene_Symbol"]+[ x for x in metap_df.columns if "CIS_" in x]
            metap_df=metap_df[gwas_columns]
            metap_df=metap_df[metap_df[f'CIS_{Disease}_LwestPvalue']<pvalue_cutof]
            metap_df=metap_df[~metap_df[f'CIS_{Disease}_LwestPvalue'].isna()]
            Exposures=list(set(metap_df[metap_df[f'CIS_{Disease}_LwestPvalue']<pvalue_cutof]["Gene_Symbol"].unique()))
            Exposures.sort()
        
        for exposure in Exposures:
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
            
            lowest_pvalue_columns=[x for x in metap_df.columns if "LwestPvalue" in x][0]
            lowest_pvalue=metap_df[metap_df["Gene_Symbol"]==exposure][lowest_pvalue_columns].apply(to_scientific).unique()[0]
            outname_final=out_name1+"_"+str(lowest_pvalue)
            
            if selected_df.shape[0]<=29:
                fp.forestplot(selected_df, estimate="b",moerror="se",ll="lo", hl="up",varlabel="SNP",decimal_precision=3,
                            ci_report=False,color_alt_rows=True,flush=False,table=True,groupvar="group",
                            rightannote=['P', 'b','up', 'lo',"nSNP","EAF"],
                            right_annoteheaders=["Pvalue", "Beta", "Up_CI","Low_CI","nSNP","EAF"],
                            sort=True,**{"marker": "D", "markersize": 35,  "xlinestyle": (0, (10, 5)), "xlinecolor": "#808080", "xtick_size": 12}   )
                plt.savefig(f"{pqtltype}_BiogenDecode_Pqtl_{Disease}/{pqtltype}_BiogenDecode_Pqtl_{Disease}_{exposure}_plot_{outname_final}.png", bbox_inches="tight")
            
            if selected_df.shape[0]>29:
                try:
                    biogen=selected_df[selected_df["SNP"].str.contains("Biogen")]
                    fp.forestplot(biogen, estimate="b",moerror="se",ll="lo", hl="up",varlabel="SNP",decimal_precision=3,
                            ci_report=False,color_alt_rows=True,flush=False,table=True,groupvar="group",
                            rightannote=['P', 'b','up', 'lo',"nSNP","EAF"],
                            right_annoteheaders=["Pvalue", "Beta", "Up_CI","Low_CI","nSNP","EAF"],
                            sort=True,**{"marker": "D", "markersize": 35,  "xlinestyle": (0, (10, 5)), "xlinecolor": "#808080", "xtick_size": 12}   )
                    plt.savefig(f"{pqtltype}_BiogenDecode_Pqtl_{Disease}/{pqtltype}_BiogenDecode_Pqtl_{Disease}_{exposure}_BiogenPqtl_plot_{outname_final}.png", bbox_inches="tight")
                except:
                    pass
                
                try:
                    decode=selected_df[selected_df["SNP"].str.contains("Decode")]
                    fp.forestplot(decode, estimate="b",moerror="se",ll="lo", hl="up",varlabel="SNP",decimal_precision=3,
                            ci_report=False,color_alt_rows=True,flush=False,table=True,groupvar="group",
                            rightannote=['P', 'b','up', 'lo',"nSNP","EAF"],
                            right_annoteheaders=["Pvalue", "Beta", "Up_CI","Low_CI","nSNP","EAF"],
                            sort=True,**{"marker": "D", "markersize": 35,  "xlinestyle": (0, (10, 5)), "xlinecolor": "#808080", "xtick_size": 12}   )
                    plt.savefig(f"{pqtltype}_BiogenDecode_Pqtl_{Disease}/{pqtltype}_BiogenDecode_Pqtl_{Disease}_{exposure}_DecodePqtl_plot_{outname_final}.png", bbox_inches="tight")
                except:
                    pass
