##This cripts for crelation plots for the MR analysis results from Decode and UKB-PPP using beta value ; Input file is Metap output file
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def plot_correlation(df, column1, column2, xlabel="X Label", ylabel="Y Label",title="Tile",disease="Diseae"):
    # Calculate the correlation between the two columns
    corr = df[column1].corr(df[column2])
    legend_properties = {'weight':'bold'}
    # Create a scatterplot
    sns.set(style="whitegrid")
    plt.figure(figsize=(8, 6))
    sns.regplot(x=column1, y=column2, data=df, ci=None, line_kws={'color':'blue'},color='red',scatter = True,scatter_kws={'s':8})  # Fit a regression line
    # Add a text label with the correlation value
    plt.legend([f'Correlation coefficient = {corr:.2f}'],loc='upper center',prop=legend_properties)
    plt.xlabel(xlabel,fontweight='bold')
    plt.ylabel(ylabel,fontweight='bold')
    plt.title(title,fontweight='bold')
    #plt.savefig(f"CIS_Correlation_plots_Decode_vs_UKB-PPP_{disease}", dpi=300, bbox_inches='tight')
    plt.savefig(f"TransExposureNoMHC_Correlation_plots_Decode_vs_UKB-PPP_{disease}", dpi=300, bbox_inches='tight')
    # Show the plot
    plt.show()


file_index="Depression_iPSYCH_2023" ## PGC3_SCZ_NoUKB ,BIP_PGC3_noukb, Cognition ,Depression_iPSYCH_2023 
disease="Depression" ## Schizophrenia "Bipolar Disorder",Cognition,Depression
decode_col = "Decode.PQTL"
biogen_col = "Biogen.PQTL"

#file=f"Biogen_Decode_pQTL_{file_index}_CisExposure_British_IVDelt_MetapAnalysis.csv"
file=f"Biogen_Decode_pQTL_{file_index}_TransExposureNoMHC_British_IVDelt_MetapAnalysis.csv"
data = pd.read_csv(file, sep=',')

column1=[x for x in data.columns if decode_col in x and 'IVWDelta_Beta' in x ][0]
column2=[x for x in data.columns if biogen_col in x and 'IVWDelta_Beta' in x ][0]

print(column2,column1)
title=f"Decode vs UKB-PPP : {disease} MR"


plot_correlation(data, column1, column2, xlabel="Beta value (Decode)", ylabel="Beta value UKB-PPP",title=title,disease=disease)
