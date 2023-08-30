import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

pd.options.display.float_format = "{:,.2f}".format
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


parser=argparse.ArgumentParser(description="This script is to create correelation altman_plot")
parser.add_argument('-inputfilename','--inputfilename', help="name of the file", required=True)
parser.add_argument('-valutype','--valutype', choices=["pvalues", "betavalues"], help="Provide the column values are P va;lue or Beta", required=True)

args = parser.parse_args()


def create_comparison_plots(data, output_bland_altman, output_correlation):
    decode_col = "Decode.PQTL"
    biogen_col = "Biogen.PQTL"
    
    decode_pvaluecol = [x for x in data.columns if decode_col in x]
    biogen_pvaluecol = [x for x in data.columns if biogen_col in x]
    
    data['NegLog_Decode_Pval'] = -np.log10(data[decode_pvaluecol[0]])
    data['NegLog_Biogen_Pval'] = -np.log10(data[biogen_pvaluecol[0]])
    data['Difference'] = data['NegLog_Decode_Pval'] - data['NegLog_Biogen_Pval']
    data['Mean'] = (data['NegLog_Decode_Pval'] + data['NegLog_Biogen_Pval']) / 2
    
    mean_diff = data['Difference'].mean()
    std_diff = data['Difference'].std()
    upper_limit = mean_diff + 1.96 * std_diff
    lower_limit = mean_diff - 1.96 * std_diff
    cr = 1.96 * std_diff
    
    plt.figure(figsize=(10, 6))
    plt.scatter(data['Mean'], data['Difference'], color='blue', s=20)
    plt.axhline(y=mean_diff, color='black', linestyle='--', label='Mean Difference')
    plt.axhline(y=upper_limit, color='red', linestyle='--', label='Upper Limit of Agreement')
    plt.axhline(y=lower_limit, color='red', linestyle='--', label='Lower Limit of Agreement')
    plt.xlabel('Mean of -log10(Decode_MR_Pipeline_pval) and -log10(Biogen_MR_Pipeline_pval)')
    plt.ylabel('Difference (-log10(Decode_MR_Pipeline_pval) - -log10(Biogen_MR_Pipeline_pval))')
    plt.text(0.05, 0.95, f'Mean Diff: {mean_diff:.4f}', transform=plt.gca().transAxes)
    plt.text(0.05, 0.90, f'Std Dev Diff: {std_diff:.4f}', transform=plt.gca().transAxes)
    plt.text(0.05, 0.85, f'Upper Limit: {upper_limit:.4f}', transform=plt.gca().transAxes)
    plt.text(0.05, 0.80, f'Lower Limit: {lower_limit:.4f}', transform=plt.gca().transAxes)
    plt.text(0.05, 0.75, f'CR: {cr:.4f}', transform=plt.gca().transAxes)
    plt.title('Bland-Altman Plot\nDifference and Limits of Agreement')
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.legend()
    plt.savefig(output_bland_altman, dpi=300, format='tiff')
    
    # Correlation Scatter Plot
    sns.set(style="whitegrid")
    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=data, x='NegLog_Decode_Pval', y='NegLog_Biogen_Pval', color='blue')
    plt.xlabel('-log10(Decode_MR_Pipeline_pval)')
    plt.ylabel('-log10(Biogen_MR_Pipeline_pval)')
    
    correlation_matrix = np.corrcoef(data['NegLog_Decode_Pval'], data['NegLog_Biogen_Pval'])
    correlation_coefficient = correlation_matrix[0, 1]
    plt.text(0.05, 0.95, f'Correlation: {correlation_coefficient:.4f}', transform=plt.gca().transAxes)
    plt.xlim([min(data['NegLog_Decode_Pval']) - 0.5, max(data['NegLog_Decode_Pval']) + 0.5])
    plt.ylim([min(data['NegLog_Biogen_Pval']) - 0.5, max(data['NegLog_Biogen_Pval']) + 0.5])
    plt.savefig(output_correlation, dpi=300, format='tiff')



def create_betabased_comparison_plots(data, output_bland_altman, output_correlation):
    decode_col = "Decode.PQTL"
    biogen_col = "Biogen.PQTL"
    
    decode_pvaluecol = [x for x in data.columns if decode_col in x]
    biogen_pvaluecol = [x for x in data.columns if biogen_col in x]
    
    data['NegLog_Decode_Pval'] = data[decode_pvaluecol[0]]
    data['NegLog_Biogen_Pval'] = data[biogen_pvaluecol[0]]
    data['Difference'] = data['NegLog_Decode_Pval'] - data['NegLog_Biogen_Pval']
    data['Mean'] = (data['NegLog_Decode_Pval'] + data['NegLog_Biogen_Pval']) / 2
    
    mean_diff = data['Difference'].mean()
    std_diff = data['Difference'].std()
    upper_limit = mean_diff + 1.96 * std_diff
    lower_limit = mean_diff - 1.96 * std_diff
    cr = 1.96 * std_diff
    
    plt.figure(figsize=(10, 6))
    plt.scatter(data['Mean'], data['Difference'], color='blue', s=20)
    plt.axhline(y=mean_diff, color='black', linestyle='--', label='Mean Difference')
    plt.axhline(y=upper_limit, color='red', linestyle='--', label='Upper Limit of Agreement')
    plt.axhline(y=lower_limit, color='red', linestyle='--', label='Lower Limit of Agreement')
    plt.xlabel('Mean of Decode Beta and Biogen Beta')
    plt.ylabel('Difference (Decode Beta - Biogen Beta)')
    plt.text(0.05, 0.95, f'Mean Diff: {mean_diff:.4f}', transform=plt.gca().transAxes)
    plt.text(0.05, 0.90, f'Std Dev Diff: {std_diff:.4f}', transform=plt.gca().transAxes)
    plt.text(0.05, 0.85, f'Upper Limit: {upper_limit:.4f}', transform=plt.gca().transAxes)
    plt.text(0.05, 0.80, f'Lower Limit: {lower_limit:.4f}', transform=plt.gca().transAxes)
    plt.text(0.05, 0.75, f'CR: {cr:.4f}', transform=plt.gca().transAxes)
    plt.title('Bland-Altman Plot\nDifference and Limits of Agreement')
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.legend()
    plt.savefig(output_bland_altman, dpi=300, format='tiff')
    
    # Correlation Scatter Plot
    sns.set(style="whitegrid")
    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=data, x='NegLog_Decode_Pval', y='NegLog_Biogen_Pval', color='blue')
    plt.xlabel('Decode Beta')
    plt.ylabel('Biogen Beta')
    
    correlation_matrix = np.corrcoef(data['NegLog_Decode_Pval'], data['NegLog_Biogen_Pval'])
    correlation_coefficient = correlation_matrix[0, 1]
    plt.text(0.05, 0.95, f'Correlation: {correlation_coefficient:.4f}', transform=plt.gca().transAxes)
    plt.xlim([min(data['NegLog_Decode_Pval']) - 0.5, max(data['NegLog_Decode_Pval']) + 0.5])
    plt.ylim([min(data['NegLog_Biogen_Pval']) - 0.5, max(data['NegLog_Biogen_Pval']) + 0.5])
    plt.savefig(output_correlation, dpi=300, format='tiff')


# Read the data from the file
if args.valutype=="pvalues":
    data = pd.read_csv(args.inputfilename, sep=',')
    create_comparison_plots(data, f'{args.inputfilename[:-4]}_altman_plot.tiff', f'{args.inputfilename[:-18]}_correlation_plot.tiff')


# Read the data from the file
if args.valutype=="betavalues":
    data = pd.read_csv(args.inputfilename, sep=',')
    create_betabased_comparison_plots(data, f'{args.inputfilename[:-4]}_altman_plot.tiff', f'{args.inputfilename[:-18]}_correlation_plot.tiff')


file="Biogen_Decode_pQTL_ASD_PGC_CisExposure_British_IVDelt_BetaValue_MetapAnalysis.csv"

data = pd.read_csv(file, sep=',')
create_betabased_comparison_plots(data, f'{file}_altman_plot.tiff', f'{file[:-18]}_correlation_plot.tiff')
