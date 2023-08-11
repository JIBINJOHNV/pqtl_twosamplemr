import argparse
import pandas as pd
import numpy as np
import os

# Create argument parser
parser = argparse.ArgumentParser(description='Process files_df')

# Add argument for files_df file path
parser.add_argument('files_df_path', type=str, help='Path to files_df file')

# Parse the command-line arguments
args = parser.parse_args()

# Read the files_df file
files_df = pd.read_csv(args.files_df_path)

# Rest of the code...
for UniqueID in files_df['UniqueID'].unique():
    result_df = pd.DataFrame(columns=['CHROM', 'GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'INFO', 'N',
                                      'TEST', 'BETA', 'SE', 'CHISQ', 'LOG10P', 'EXTRA', 'MAF', 'GeneSymbol',
                                      'UNIPROT_ID', 'CAtegary', 'GENPOS_Hg19', 'rsid'])
    gene_df = files_df[files_df['UniqueID'] == UniqueID].reset_index().drop("index", axis=1)
    gene_df=gene_df.sort_values(by="file_name")
    for Index in range(0, gene_df.shape[0]):
        df = pd.read_csv(gene_df.iloc[Index, 0], sep="\s")
        df['MAF'] = np.where(df['A1FREQ'] < 0.5, df['A1FREQ'], 1 - df['A1FREQ'])
        df['GeneSymbol'] = gene_df.loc[Index, "GeneSymbol"]
        df['UNIPROT_ID'] = gene_df.loc[Index, "UNIPROT_ID"]
        df['CAtegary'] = gene_df.loc[Index, "CAtegary"]
        df['GENPOS_Hg19'] = df['ID'].str.split(":", expand=True)[1]
        df["rsid"] = df['CHROM'].astype(str) + ":" + df['GENPOS_Hg19'].astype(str)
        result_df = pd.concat([result_df, df])
    result_df.to_csv(f'discovery_RS2_{gene_df.loc[Index, "GeneSymbol"]}:{gene_df.loc[Index, "UNIPROT_ID"]}:{gene_df.loc[Index, "CAtegary"]}.tsv', sep="\t")
    os.system(f'gzip discovery_RS2_{gene_df.loc[Index, "GeneSymbol"]}:{gene_df.loc[Index, "UNIPROT_ID"]}:{gene_df.loc[Index, "CAtegary"]}.tsv')
    print(result_df)


#GWAS summary statistics located in /edgehpc/dept/compbio/human_genetics/users/bsun/UKB_PPP_pQTLs/analysis_v2/interim/UKBPPP_GWAS_v2/RS2/
