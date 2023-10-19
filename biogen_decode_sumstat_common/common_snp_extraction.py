
import pandas as pd
from scipy.stats import pearsonr


decode="/edgehpc/dept/human_genetics/pQTL_deCODE_Ferkingstad_2021/gwas/7933_75_ADAM22_ADA22.txt.gz"
biogen="/edgehpc/dept/human_genetics/users/jjohn1/Discovery_Genewise_Merged/discovery_RS2_ADAM22:Q9P0K1:Neurology.tsv.gz"

cor_dict={}
sig_cor_dict={}
sig_betacor_dict={}
maf01_sig_cor_dict={}
maf01_sig_betacor_dict={}


d_df=pd.read_csv(decode,sep="\t")
b_df=pd.read_csv(biogen,sep="\t")

b_df["CHROM"]="chr"+b_df["CHROM"].astype("str")

b_df2=b_df[['CHROM','GENPOS','ALLELE0', 'ALLELE1','BETA', 'SE','LOG10P','CAtegary','A1FREQ',"MAF","INFO","N","GeneSymbol"]]
d_df2=d_df[['Chrom', 'Pos','effectAllele', 'otherAllele','Beta','SE','minus_log10_pval',"N","ImpMAF"]]
b_df2['GENPOS']=b_df2['GENPOS'].astype("int")
d_df2['Pos']=d_df2['Pos'].astype("int")

b_df2=b_df2[ (b_df2["ALLELE0"].str.len()==1) & (b_df2["ALLELE1"].str.len()==1) & (b_df2["MAF"]>0.001) ]
d_df2=d_df2[ (d_df2["effectAllele"].str.len()==1) & (d_df2["otherAllele"].str.len()==1) & (d_df2["ImpMAF"]>0.001) ]


##Add prefix to column names
b_df2.columns=["biogen_"+x for x in b_df2.columns]
d_df2.columns=["decode_"+x for x in d_df2.columns]


bi_de_df1=pd.merge(b_df2,d_df2,left_on=["biogen_CHROM","biogen_GENPOS","biogen_ALLELE1","biogen_ALLELE0"],right_on=["decode_Chrom","decode_Pos",'decode_effectAllele', 'decode_otherAllele'])
bi_de_df2=pd.merge(b_df2,d_df2,left_on=["biogen_CHROM","biogen_GENPOS","biogen_ALLELE1","biogen_ALLELE0"],right_on=["decode_Chrom","decode_Pos",'decode_otherAllele','decode_effectAllele', ])
bi_de_df3=pd.merge(b_df2,d_df2,left_on=["biogen_CHROM","biogen_GENPOS","biogen_ALLELE0","biogen_ALLELE1"],right_on=["decode_Chrom","decode_Pos",'decode_effectAllele', 'decode_otherAllele'])


bi_de_df=pd.concat([bi_de_df1,bi_de_df2,bi_de_df3])

correlation_coefficient, p_value = pearsonr(bi_de_df['biogen_LOG10P'],bi_de_df['decode_minus_log10_pval'])
gene=bi_de_df.head(n=1)['biogen_GeneSymbol'].iloc[0]
cor_dict[gene]=[correlation_coefficient, p_value]


biogen_columns=[x for x in bi_de_df.columns if "biogen_" in x]
decode_columns=[x for x in bi_de_df.columns if "decode_" in x]+['biogen_GeneSymbol']

#biogen_df=bi_de_df[biogen_columns].rename(columns={"biogen_GeneSymbol":"GeneSymbol"})
#decode_df=bi_de_df[decode_columns].rename(columns={"biogen_GeneSymbol":"GeneSymbol"})

sig_001= bi_de_df[(bi_de_df['biogen_LOG10P']>5) | (bi_de_df['decode_minus_log10_pval']>5)]

correlation_coefficient, p_value = pearsonr(sig_001['biogen_LOG10P'],sig_001['decode_minus_log10_pval'])
sig_cor_dict[gene]=[correlation_coefficient, p_value]

correlation_coefficient, p_value = pearsonr(sig_001['biogen_BETA'],sig_001['decode_Beta'])
sig_betacor_dict[gene]=[correlation_coefficient, p_value]

biogen_p001df=sig_001[biogen_columns].rename(columns={"biogen_GeneSymbol":"GeneSymbol"})
decode_p001df=sig_001[decode_columns].rename(columns={"biogen_GeneSymbol":"GeneSymbol"})


biogen_p001df.to_csv("biogen_gene_MAF001.tsv",sep="\t")
decode_p001df.to_csv("biogen_gene_MAF001.tsv",sep="\t")



sig_01=sig_001[(sig_001['biogen_MAF']>0.01) & (sig_001['decode_ImpMAF']>0.01) ]

correlation_coefficient, p_value = pearsonr(sig_01['biogen_LOG10P'],sig_01['decode_minus_log10_pval'])
maf01_sig_cor_dict[gene]=[correlation_coefficient, p_value]

correlation_coefficient, p_value = pearsonr(sig_01['biogen_BETA'],sig_01['decode_Beta'])
maf01_sig_betacor_dict[gene]=[correlation_coefficient, p_value]


biogen_p01df=sig_01[biogen_columns].rename(columns={"biogen_GeneSymbol":"GeneSymbol"})
decode_p01df=sig_01[decode_columns].rename(columns={"biogen_GeneSymbol":"GeneSymbol"})


biogen_p01df.to_csv("biogen_gene_MAF001.tsv",sep="\t")
decode_p01df.to_csv("biogen_gene_MAF001.tsv",sep="\t")


