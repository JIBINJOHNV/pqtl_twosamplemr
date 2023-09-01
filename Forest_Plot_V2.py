import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


biogen_excell="Biogen_ASD_PGC_TransExposureNoMHC_CompleteMR_AnalysisResults.xlsx"
decode_exsell="Decode_ASD_PGC_TransExposureNoMHC_CompleteMR_AnalysisResults.xlsx"
metap_asd="TransExposureNoMHCExposure_CompleteMR.xlsx"
Disease="ASD_PGC"

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
    
    return mr_df2


def  single_mr(pqtl_type,excell):
    mr_df=pd.read_excel(excell,sheet_name="Singlevariant")
    selected_columns=['Gene_Symbol','outcome','TwosampleMR-wald_Beta', 'TwosampleMR-wald_SE','TwosampleMR-wald_Pvalue',"SNP"]
    mr_df2=mr_df[selected_columns]
    mr_df2.columns=['id.exposure', 'id.outcome','b', 'se','P',"SNP"]
    mr_df2['up'] = mr_df2['b'] + 1.96 * mr_df2['se']
    mr_df2['lo'] = mr_df2['b'] - 1.96 * mr_df2['se']
    
    if pqtl_type=="biogen":
        mr_df2["group"]="Biogen_SinglSNP"
        mr_df2["SNP"]="Biogen_"+mr_df2["SNP"]
    elif pqtl_type=="decode":
        mr_df2["group"]="Decode_SinglSNP"
        mr_df2["SNP"]="Decode_"+mr_df2["SNP"]
    
    mr_df2=mr_df2.sort_values(by=['id.exposure',"SNP"])
    mr_df2['nSNP']="_"
    return mr_df2



biogen_mr_df=combined_mr("biogen",biogen_excell)
decode_mr_df=combined_mr("decode",decode_exsell)

biogen_single_df=single_mr("biogen",biogen_excell)
decode_single_df=single_mr("decode",decode_exsell)

mr_df=pd.concat([biogen_mr_df,decode_mr_df])
mr_single=pd.concat([biogen_single_df,decode_single_df])



metap_df=pd.read_excel(metap_asd)
gwas_columns=["Gene_Symbol"]+[ x for x in metap_df.columns if Disease in x]
metap_df=metap_df[gwas_columns]
metap_df=metap_df[metap_df[f'TransExposureNoMHC_{Disease}_LwestPvalue']<pvalue_cutof]
metap_df=metap_df[~metap_df[f'TransExposureNoMHC_{Disease}_LwestPvalue'].isna()]
Exposures=list(metap_df[metap_df[f'TransExposureNoMHC_{Disease}_LwestPvalue']<pvalue_cutof]["Gene_Symbol"].unique())
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
    
    if selected_df.shape[0]<=29:
        fp.forestplot(selected_df, estimate="b",moerror="se",ll="lo", hl="up",varlabel="SNP",decimal_precision=3,
                    ci_report=False,color_alt_rows=True,flush=False,table=True,groupvar="group",
                    rightannote=['P', 'b','up', 'lo',"nSNP"],
                    right_annoteheaders=["Pvalue", "Beta", "Up_CI","Low_CI","nSNP"],
                    sort=True,**{"marker": "D", "markersize": 35,  "xlinestyle": (0, (10, 5)), "xlinecolor": "#808080", "xtick_size": 12}   )
        plt.savefig(f"TransNoMHC_BiogenDecode_Pqtl_{Disease}_{exposure}_plot.png", bbox_inches="tight")
    
    if selected_df.shape[0]>=29:
        try:
            biogen=selected_df[selected_df["SNP"].str.contains("Biogen")]
            
            fp.forestplot(biogen, estimate="b",moerror="se",ll="lo", hl="up",varlabel="SNP",decimal_precision=3,
                    ci_report=False,color_alt_rows=True,flush=False,table=True,groupvar="group",
                    rightannote=['P', 'b','up', 'lo',"nSNP"],
                    right_annoteheaders=["Pvalue", "Beta", "Up_CI","Low_CI","nSNP"],
                    sort=True,**{"marker": "D", "markersize": 35,  "xlinestyle": (0, (10, 5)), "xlinecolor": "#808080", "xtick_size": 12}   )
            plt.savefig(f"TransNoMHC_BiogenDecode_Pqtl_{Disease}_{exposure}_BiogenPqtl_plot.png", bbox_inches="tight")
        except:
            pass
        
        try:
            decode=selected_df[selected_df["SNP"].str.contains("Decode")]
            fp.forestplot(decode, estimate="b",moerror="se",ll="lo", hl="up",varlabel="SNP",decimal_precision=3,
                    ci_report=False,color_alt_rows=True,flush=False,table=True,groupvar="group",
                    rightannote=['P', 'b','up', 'lo',"nSNP"],
                    right_annoteheaders=["Pvalue", "Beta", "Up_CI","Low_CI","nSNP"],
                    sort=True,**{"marker": "D", "markersize": 35,  "xlinestyle": (0, (10, 5)), "xlinecolor": "#808080", "xtick_size": 12}   )
            plt.savefig(f"TransNoMHC_BiogenDecode_Pqtl_{Disease}_{exposure}_DecodePqtl_plot.png", bbox_inches="tight")
        except:
            pass


