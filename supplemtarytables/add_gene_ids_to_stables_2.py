import numpy as np
import pandas as pd


mr_decode=pd.read_csv("Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2.csv")
mr_decode=mr_decode[['Protein','DECODE_exposure']]
mr_decode=mr_decode[mr_decode["DECODE_exposure"].notna()].drop_duplicates()
mr_decode["seqid.1"]=mr_decode["DECODE_exposure"].str.split("_",expand=True)[1]+"_"+mr_decode["DECODE_exposure"].str.split("_",expand=True)[2]
decod=pd.read_excel("olink_somscan_decode_biobank_protein.xlsx",sheet_name="Sheet1")
full_decode=pd.merge(mr_decode,decod,on="seqid.1",how="left")


##
biomart=pd.read_csv("biomart_gene_details_109.csv")
biomart[['HGNC','NCBI_Gene_ID','OMIM_Gene_ID']]=biomart[['HGNC','NCBI_Gene_ID','OMIM_Gene_ID']].astype('Int64')



###### HGNC
hgnc=pd.read_csv("hgnc_custom.tsv",sep="\t")
hgnc["HGNC ID"]=hgnc["HGNC ID"].str.split(":",expand=True)[1]
hgnc=hgnc.rename(columns={"HGNC ID":"HGNC","NCBI Gene ID":"GeneID_hgnc"})
hgnc[["HGNC","GeneID_hgnc","NCBI Gene ID(supplied by NCBI)"]]=hgnc[["HGNC","GeneID_hgnc","NCBI Gene ID(supplied by NCBI)"]].astype('Int64')
hgnc["Ensembl_gene_ID"]=np.where(hgnc["Ensembl gene ID"].isna(),hgnc["Ensembl gene ID"],hgnc["Ensembl ID(supplied by Ensembl)"])
hgnc["NCBI_Gene_ID"]=np.where(hgnc["GeneID_hgnc"].isna(),hgnc["NCBI Gene ID(supplied by NCBI)"],hgnc["GeneID_hgnc"])
hgnc=hgnc.drop(["GeneID_hgnc","NCBI Gene ID(supplied by NCBI)","Ensembl gene ID","Ensembl ID(supplied by Ensembl)","Accession numbers"],axis=1)
hgnc=hgnc.rename(columns={'UniProt ID(supplied by UniProt)':"UniProt_ID",'OMIM ID(supplied by OMIM)':'OMIM_Gene_ID'})
hgnc[['Previous symbols', 'Alias symbols', 'Approved symbol']] = hgnc[['Previous symbols', 'Alias symbols', 'Approved symbol']].apply(lambda x: x.str.upper())

hgn_previousysmbols=hgnc.copy()
hgn_previousysmbols['Previous symbols']=hgn_previousysmbols['Previous symbols'].str.split(",")
hgn_previousysmbols=hgn_previousysmbols.explode('Previous symbols')
hgn_previousysmbols['Previous symbols']=hgn_previousysmbols['Previous symbols'].str.strip()

hgnc_alias=hgnc.copy()
hgnc_alias['Alias symbols']=hgnc_alias['Alias symbols'].str.split(",")
hgnc_alias=hgnc_alias.explode('Alias symbols')
hgnc_alias['Alias symbols']=hgnc_alias['Alias symbols'].str.strip()


##Merging with approveedSymbols
## Approved symbols
full_decode_hgnc1=pd.merge(full_decode,hgnc,left_on=["uniprot","gene_name"],right_on=["UniProt_ID","Approved symbol"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_approveed_genename=full_decode_hgnc1[full_decode_hgnc1["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
approveed_genename=full_decode_hgnc1[~full_decode_hgnc1["Approved symbol"].isna()]

full_decode_hgnc2=pd.merge(no_approveed_genename,hgnc,left_on=["uniprot","target_name"],right_on=["UniProt_ID","Approved symbol"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_approveed_tname=full_decode_hgnc2[full_decode_hgnc2["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
approveed_tname=full_decode_hgnc2[~full_decode_hgnc2["Approved symbol"].isna()]

##Previous symbols
full_decode_hgnc3=pd.merge(no_approveed_tname,hgn_previousysmbols,left_on=["uniprot","gene_name"],right_on=["UniProt_ID","Previous symbols"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_previous_genename=full_decode_hgnc3[full_decode_hgnc3["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
previous_genename=full_decode_hgnc3[~full_decode_hgnc3["Approved symbol"].isna()]

full_decode_hgnc4=pd.merge(no_previous_genename,hgn_previousysmbols,left_on=["uniprot","target_name"],right_on=["UniProt_ID","Previous symbols"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_previous_tname=full_decode_hgnc4[full_decode_hgnc4["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
previous_tname=full_decode_hgnc4[~full_decode_hgnc4["Approved symbol"].isna()]

##Alias symbols
full_decode_hgnc5=pd.merge(no_previous_tname,hgnc_alias,left_on=["uniprot","gene_name"],right_on=["UniProt_ID","Alias symbols"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_alias_genename=full_decode_hgnc5[full_decode_hgnc5["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
alias_genename=full_decode_hgnc5[~full_decode_hgnc5["Approved symbol"].isna()]

full_decode_hgnc4=pd.merge(no_alias_genename,hgnc_alias,left_on=["uniprot","target_name"],right_on=["UniProt_ID","Alias symbols"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_alias_tname=full_decode_hgnc5[full_decode_hgnc5["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
alias_tname=full_decode_hgnc5[~full_decode_hgnc5["Approved symbol"].isna()]


########Biomart
full_decode_biomart1=pd.merge(no_alias_tname,biomart,left_on=["uniprot","Protein"],right_on=["UniProt_ID","Approved symbol"],how="left").drop(["Gene Synonym",'Gene name'],axis=1).drop_duplicates()
no_biomart_approvedsymobols=full_decode_biomart1[full_decode_biomart1["Approved symbol"].isna()].drop(['HGNC',"Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
biomart_approvedsymobols=full_decode_biomart1[~full_decode_biomart1["Approved symbol"].isna()]

##remaining from uniprot
not_hgnc_biomart=[{"UniProt_ID":"P04745","Approved symbol":"AMY1A","Ensembl_gene_ID":"ENSG00000237763","NCBI_Gene_ID":276,"HGNC":474,"OMIM_Gene_ID":104700},
 {"UniProt_ID":"J3QR46","Approved symbol":"KIAA0040","Ensembl_gene_ID":"ENSG00000235750","NCBI_Gene_ID":9674,"HGNC":28950,"OMIM_Gene_ID":616696},
  {"UniProt_ID":"Q9Y4X1","Approved symbol":"UGT2A1","Ensembl_gene_ID":"ENSG00000173610","NCBI_Gene_ID":10941,"HGNC":12542,"OMIM_Gene_ID":604716},
  {"UniProt_ID":"Q8NFS9","Approved symbol":"GCNT2","Ensembl_gene_ID":"ENSG00000111846","NCBI_Gene_ID":2651,"HGNC":4204,"OMIM_Gene_ID":600429}]
not_hgnc_biomart_df=pd.DataFrame(not_hgnc_biomart)
full_decode_custom=pd.merge(no_biomart_approvedsymobols,not_hgnc_biomart_df,left_on="uniprot",right_on="UniProt_ID")


decode_withgeneinfo=pd.concat([approveed_genename,approveed_tname,previous_genename,previous_tname,
                               alias_genename,alias_tname,biomart_approvedsymobols,full_decode_custom ]).drop_duplicates().drop("Approved name",axis=1)

decode_withgeneinfo=decode_withgeneinfo[['Protein','DECODE_exposure','seqid.1', 'seqid','platform', 'seqid2','platform','target_full_name',
                        'target_name','gene_name','uniprot','UniProt_ID','Approved symbol','HGNC','OMIM_Gene_ID','Ensembl_gene_ID', 'NCBI_Gene_ID']]

decode_withgeneinfo.to_csv("Complete_Decodee_protein_withcomplete_Geneinformation.csv",index=None)
############------------------------------------------############################################################################

ukp_2=pd.read_csv("ukbb-ppp_protein_details.txt",sep=" ",header=None)
ukp_2.columns=["gene_name","uniprot","oid","version","panel"]
ukp_1=pd.read_excel("olink_somscan_decode_biobank_protein.xlsx",sheet_name="Sheet2")
ukpp=pd.merge(ukp_1,ukp_2,on=["gene_name","uniprot","oid"])

mr_decode=pd.read_csv("Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table2.csv")[["Protein","UKBB_PPP_Panel",'UKBB_PPP_exposure']].drop_duplicates()
mr_decode=mr_decode[~mr_decode["UKBB_PPP_exposure"].isna()].drop_duplicates()

full_ukbb_ppp=pd.merge(mr_decode,ukpp,left_on="Protein",right_on="gene_name",how="left").drop_duplicates()
full_ukbb_ppp["panel"]=full_ukbb_ppp["panel"]+"_"+full_ukbb_ppp["version"]
full_ukbb_ppp=full_ukbb_ppp.drop("version",axis=1)
full_ukbb_ppp.to_csv("Complete_ukbb_ppp_protein_olink.csv",index=None)

full_ukbb_ppp=pd.read_csv("Complete_ukbb_ppp_protein_olink.csv")



##Merging with approveedSymbols
## Approved symbols
full_decode=full_ukbb_ppp.copy()
full_decode=full_decode.fillna("NA")
full_decode_hgnc1=pd.merge(full_decode,hgnc,left_on=["uniprot","gene_name"],right_on=["UniProt_ID","Approved symbol"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_approveed_genename=full_decode_hgnc1[full_decode_hgnc1["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
approveed_genename=full_decode_hgnc1[~full_decode_hgnc1["Approved symbol"].isna()]

full_decode_hgnc2=pd.merge(no_approveed_genename,hgnc,left_on=["uniprot","Protein"],right_on=["UniProt_ID","Approved symbol"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_approveed_tname=full_decode_hgnc2[full_decode_hgnc2["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
approveed_tname=full_decode_hgnc2[~full_decode_hgnc2["Approved symbol"].isna()]

##Previous symbols
full_decode_hgnc3=pd.merge(no_approveed_tname,hgn_previousysmbols,left_on=["uniprot","gene_name"],right_on=["UniProt_ID","Previous symbols"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_previous_genename=full_decode_hgnc3[full_decode_hgnc3["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
previous_genename=full_decode_hgnc3[~full_decode_hgnc3["Approved symbol"].isna()]

full_decode_hgnc4=pd.merge(no_previous_genename,hgn_previousysmbols,left_on=["uniprot","Protein"],right_on=["UniProt_ID","Previous symbols"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_previous_tname=full_decode_hgnc4[full_decode_hgnc4["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
previous_tname=full_decode_hgnc4[~full_decode_hgnc4["Approved symbol"].isna()]

##Alias symbols
full_decode_hgnc5=pd.merge(no_previous_tname,hgnc_alias,left_on=["uniprot","gene_name"],right_on=["UniProt_ID","Alias symbols"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_alias_genename=full_decode_hgnc5[full_decode_hgnc5["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
alias_genename=full_decode_hgnc5[~full_decode_hgnc5["Approved symbol"].isna()]

full_decode_hgnc6=pd.merge(no_alias_genename,hgnc_alias,left_on=["uniprot","Protein"],right_on=["UniProt_ID","Alias symbols"],how="left").drop(["Previous symbols","Alias symbols","Status"],axis=1).drop_duplicates()
no_alias_tname=full_decode_hgnc6[full_decode_hgnc6["Approved symbol"].isna()].drop(['HGNC','Approved name', "Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
alias_tname=full_decode_hgnc6[~full_decode_hgnc6["Approved symbol"].isna()]



########Biomart
full_decode_biomart1=pd.merge(no_alias_tname,biomart,left_on=["uniprot","Protein"],right_on=["UniProt_ID","Approved symbol"],how="left").drop(["Gene Synonym",'Gene name'],axis=1).drop_duplicates()
no_biomart_approvedsymobols=full_decode_biomart1[full_decode_biomart1["Approved symbol"].isna()].drop(['HGNC',"Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
biomart_approvedsymobols=full_decode_biomart1[~full_decode_biomart1["Approved symbol"].isna()]


biomart2=biomart[biomart["Approved symbol"].isin(["PIBF1","ATRN","TRAF3","BCL2L11","NTproBNP"])].drop("Gene Synonym",axis=1).drop_duplicates()
biomart2=biomart2[~biomart2["UniProt_ID"].isna()]

full_decode_biomart2=pd.merge(no_alias_tname,biomart2,left_on=["Protein"],right_on=["Approved symbol"],how="left").drop(['Gene name'],axis=1).drop_duplicates()
no_biomart_approvedsymobols2=full_decode_biomart2[full_decode_biomart2["Approved symbol"].isna()].drop(['HGNC',"Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
biomart_approvedsymobols2=full_decode_biomart2[~full_decode_biomart2["Approved symbol"].isna()]


full_decode_biomart3=pd.merge(no_biomart_approvedsymobols2,biomart,left_on=["uniprot"],right_on=["UniProt_ID"],how="left").drop(['Gene name',"Gene Synonym"],axis=1).drop_duplicates()
no_biomart_approvedsymobols3=full_decode_biomart3[full_decode_biomart3["Approved symbol"].isna()].drop(['HGNC',"Approved symbol",'UniProt_ID','OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID'],axis=1)
biomart_approvedsymobols3=full_decode_biomart3[~full_decode_biomart3["Approved symbol"].isna()]
biomart_approvedsymobols3=biomart_approvedsymobols3.groupby(['Protein', 'UKBB_PPP_Panel', 'UKBB_PPP_exposure', 'oid', 'uniprot', 'target_full_name', 
                                        'gene_name', 'platform', 'panel','HGNC', 'Approved symbol', 'NCBI_Gene_ID', 'UniProt_ID', 'OMIM_Gene_ID'])['Ensembl_gene_ID'].agg(lambda x: ', '.join(x)).reset_index()





ukbb_ppp_withgeneinfo=pd.concat([approveed_genename,approveed_tname,previous_genename,previous_tname,alias_genename,alias_tname,biomart_approvedsymobols,
                                 biomart_approvedsymobols2,biomart_approvedsymobols3,no_biomart_approvedsymobols3]).drop_duplicates()
ukbb_ppp_withgeneinfo=ukbb_ppp_withgeneinfo[~ukbb_ppp_withgeneinfo['Ensembl_gene_ID'].isin(['ENSG00000274287','ENSG00000180900'])]

ukbb_ppp_withgeneinfo=ukbb_ppp_withgeneinfo.rename(columns={'oid':'UKBB_PPP_oid','target_full_name':'UKBB_PPP_target_full_name','gene_name':'UKBB_PPP_gene_name',
                                      'platform':'UKBB_PPP_platform','panel':'UKBB_PPP_panel','uniprot':'UKBB_PPP_uniprot','Protein':'UKBB_PPP_Protein'})
ukbb_ppp_withgeneinfo=ukbb_ppp_withgeneinfo.drop('Approved name',axis=1)
ukbb_ppp_withgeneinfo[[x for x in ukbb_ppp_withgeneinfo.columns if "HGNC" not in  x]]=ukbb_ppp_withgeneinfo[[x for x in ukbb_ppp_withgeneinfo.columns if "HGNC" not in  x]].fillna("NA")

ukbb_ppp_withgeneinfo.to_csv("Complete_UKBB_PPPP_protein_withcomplete_Geneinformation.csv",index=None)


decode_withgeneinfo=decode_withgeneinfo.rename(columns={'seqid.1':'DECODE_seqid.1','seqid':'DECODE_seqid','platform':'DECODE_platform','seqid2':'DECODE_seqid2',
                          'platform':'DECODE_platform','target_full_name':'DECODE_target_full_name','target_name':'DECODE_target_name',
                          'gene_name':'DECODE_gene_name','uniprot':'DECODE_uniprot','Protein':'DECODE_Protein'})
columns_to_fill = [col for col in decode_withgeneinfo.columns if "HGNC" not in col]
decode_withgeneinfo[columns_to_fill].fillna("NA",inplace=True)

decode_withgeneinfo.to_csv("Complete_Decodee_protein_withcomplete_Geneinformation.csv",index=None)


decode_withgeneinfo['OMIM_Gene_ID']=decode_withgeneinfo['OMIM_Gene_ID'].astype("str".strip())
ukbb_ppp_withgeneinfo['OMIM_Gene_ID']=ukbb_ppp_withgeneinfo['OMIM_Gene_ID'].astype("str".strip())

ukbb_decode_withgeneinfo=pd.merge(decode_withgeneinfo,ukbb_ppp_withgeneinfo,on=["HGNC","Approved symbol","UniProt_ID",'NCBI_Gene_ID','Ensembl_gene_ID'],how="outer")
ukbb_decode_withgeneinfo=ukbb_decode_withgeneinfo.drop(['DECODE_seqid2','DECODE_seqid.1'],axis=1)
ukbb_decode_withgeneinfo.to_csv("Complete_Decodee_UKBB_PPP_protein_withcomplete_Geneinformation.csv",index=None)
