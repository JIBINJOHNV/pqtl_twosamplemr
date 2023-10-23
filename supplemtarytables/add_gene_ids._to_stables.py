##################---------------------------------------Gene ID's to UKBB-PPP and Deecode --------------- ##########################################
decod=pd.read_excel("olink_somscan_decode_biobank_protein.xlsx",sheet_name="Sheet1")
ukpp=pd.read_excel("olink_somscan_decode_biobank_protein.xlsx",sheet_name="Sheet2")

mr=pd.read_excel("olink_somscan_decode_biobank_protein.xlsx",sheet_name="Sheet3")
mr_decode=mr[~mr["DECODE_exposure"].isna()][["Protein","DECODE_exposure"]].drop_duplicates()
mr_ukbppp=mr[~mr["UKBB_PPP_exposure"].isna()][["Protein","UKBB_PPP_exposure"]].drop_duplicates()

mr_decode_decodeoriginal=pd.merge(mr_decode,decod,left_on=["Protein","DECODE_exposure"],right_on=["gene_name","seqid2"],how="left")
mr_decode_decodeoriginal.to_csv("mr_decode_decodeoriginal.csv",index=None)


###### HGNC
hgnc=pd.read_csv("hgnc_custom.tsv",sep="\t")
hgnc["HGNC ID"]=hgnc["HGNC ID"].str.split(":",expand=True)[1]
hgnc=hgnc.rename(columns={"HGNC ID":"HGNC","NCBI Gene ID":"GeneID_hgnc"})
hgnc[["HGNC","GeneID_hgnc","NCBI Gene ID(supplied by NCBI)"]]=hgnc[["HGNC","GeneID_hgnc","NCBI Gene ID(supplied by NCBI)"]].astype('Int64')
hgnc["Ensembl_gene_ID"]=np.where(hgnc["Ensembl gene ID"].isna(),hgnc["Ensembl gene ID"],hgnc["Ensembl ID(supplied by Ensembl)"])
hgnc["NCBI_Gene_ID"]=np.where(hgnc["GeneID_hgnc"].isna(),hgnc["NCBI Gene ID(supplied by NCBI)"],hgnc["GeneID_hgnc"])
hgnc=hgnc.drop(["GeneID_hgnc","NCBI Gene ID(supplied by NCBI)","Ensembl gene ID","Ensembl ID(supplied by Ensembl)","Accession numbers"],axis=1)
hgnc=hgnc.rename(columns={'UniProt ID(supplied by UniProt)':"UniProt_ID",'OMIM ID(supplied by OMIM)':'OMIM_Gene_ID'})
#hgnc['Alias symbols']=hgnc['Alias symbols'].str.replace(" ", "")
#hgnc['Previous symbols']=hgnc['Previous symbols'].str.replace(" ", "")

hgn_previousysmbols=hgnc.copy()
hgn_previousysmbols['Previous symbols']=hgn_previousysmbols['Previous symbols'].str.split(",")
hgn_previousysmbols=hgn_previousysmbols.explode('Previous symbols')
hgn_previousysmbols['Previous symbols']=hgn_previousysmbols['Previous symbols'].str.strip()

hgnc_alias=hgnc.copy()
hgnc_alias['Alias symbols']=hgnc_alias['Alias symbols'].str.split(",")
hgnc_alias=hgnc_alias.explode('Alias symbols')
hgnc_alias['Alias symbols']=hgnc_alias['Alias symbols'].str.strip()

mr_decode_decodeoriginal_hgnc=pd.merge(mr_decode_decodeoriginal,hgnc,left_on=["uniprot","Protein"],right_on=["UniProt_ID","Approved symbol"],how="left").drop_duplicates()
not_approved=mr_decode_decodeoriginal_hgnc[mr_decode_decodeoriginal_hgnc["UniProt_ID"].isna()].iloc[:,0:10]
approved_symbols=mr_decode_decodeoriginal_hgnc[mr_decode_decodeoriginal_hgnc["UniProt_ID"].notna()]

mr_decode_decodeoriginal_hgnc_previous=pd.merge(not_approved,hgn_previousysmbols,left_on=["uniprot","Protein"],right_on=["UniProt_ID","Previous symbols"],how="left").drop_duplicates()
not_previous=mr_decode_decodeoriginal_hgnc_previous[mr_decode_decodeoriginal_hgnc_previous["UniProt_ID"].isna()].iloc[:,0:10]
previous_symbols=mr_decode_decodeoriginal_hgnc_previous[mr_decode_decodeoriginal_hgnc_previous["UniProt_ID"].notna()]

mr_decode_decodeoriginal_hgnc_alias=pd.merge(not_previous,hgnc_alias,left_on=["uniprot","Protein"],right_on=["UniProt_ID",'Alias symbols'],how="left").drop_duplicates()
not_alias=mr_decode_decodeoriginal_hgnc_alias[mr_decode_decodeoriginal_hgnc_alias["UniProt_ID"].isna()].iloc[:,0:10]
alias_symbols=mr_decode_decodeoriginal_hgnc_previous[mr_decode_decodeoriginal_hgnc_previous["UniProt_ID"].notna()]

biomart=pd.read_csv("/Users/jibinjohn/Desktop/MR_Analysis_supplementarytablees/biomart_gene_details_109.csv")
biomart[['HGNC','NCBI_Gene_ID','OMIM_Gene_ID']]=biomart[['HGNC','NCBI_Gene_ID','OMIM_Gene_ID']].astype('Int64')
biomart_1=biomart.drop(["Gene Synonym",'Gene name'],axis=1).drop_duplicates()
biomart_2=biomart.drop(['Gene name'],axis=1).drop_duplicates()


mr_decode_decodeoriginal_biomart_approved=pd.merge(not_alias,biomart_1,left_on=["uniprot","Protein"],right_on=["UniProt_ID",'Approved symbol'],how="left").drop_duplicates()
not_biomart_approved=mr_decode_decodeoriginal_biomart_approved[mr_decode_decodeoriginal_biomart_approved["Approved symbol"].isna()].iloc[:,0:10]
biomart_approved=mr_decode_decodeoriginal_biomart_approved[mr_decode_decodeoriginal_biomart_approved["Approved symbol"].notna()]

mr_decode_decodeoriginal_biomart_synonym=pd.merge(not_biomart_approved,biomart_2,left_on=["uniprot","Protein"],right_on=["UniProt_ID",'Gene Synonym'],how="left").drop_duplicates()
not_biomart_synonym=mr_decode_decodeoriginal_biomart_synonym[mr_decode_decodeoriginal_biomart_synonym["Approved symbol"].isna()].iloc[:,0:10]
biomart_synonym=mr_decode_decodeoriginal_biomart_synonym[mr_decode_decodeoriginal_biomart_synonym["Approved symbol"].notna()]
biomart_synonym=biomart_synonym[biomart_synonym["Ensembl_gene_ID"].isin(["ENSG00000112667","ENSG00000166503",
                                                    "ENSG00000204420","ENSG00000167674"])].drop_duplicates().sort_values(by="gene_name")

not_hgnc_biomart=[{"UniProt_ID":"P04745","Approved symbol":"AMY1A","Ensembl_gene_ID":"ENSG00000237763","NCBI_Gene_ID":276,"HGNC":474,"OMIM_Gene_ID":104700},
 {"UniProt_ID":"J3QR46","Approved symbol":"KIAA0040","Ensembl_gene_ID":"ENSG00000235750","NCBI_Gene_ID":9674,"HGNC":28950,"OMIM_Gene_ID":616696},
  {"UniProt_ID":"Q9Y4X1","Approved symbol":"UGT2A1","Ensembl_gene_ID":"ENSG00000173610","NCBI_Gene_ID":10941,"HGNC":12542,"OMIM_Gene_ID":604716},
  {"UniProt_ID":"Q8NFS9","Approved symbol":"GCNT2","Ensembl_gene_ID":"ENSG00000111846","NCBI_Gene_ID":2651,"HGNC":4204,"OMIM_Gene_ID":600429}]

not_hgnc_biomart_df=pd.DataFrame(not_hgnc_biomart)
not_biomart_synonym_custom_df=pd.merge(not_biomart_synonym,not_hgnc_biomart_df,left_on="uniprot",right_on="UniProt_ID")

decode_full=pd.concat([approved_symbols,previous_symbols,alias_symbols,biomart_approved,
                                        biomart_synonym,not_biomart_synonym_custom_df]).drop(["Alias symbols","Gene Synonym","Previous symbols"],axis=1).drop_duplicates()


decode_full.to_csv("Complete_Decodee_protein_withcomplete_Geneinformation.csv",index=None)

#######################---------------------------UKBB-PPP-------------------------------######################################

ukp_2=pd.read_csv("ukbb-ppp_protein_details.txt",sep=" ",header=None)
ukp_2.columns=["gene_name","uniprot","oid","version","panel"]
ukp_1=pd.read_excel("olink_somscan_decode_biobank_protein.xlsx",sheet_name="Sheet2")
ukpp=pd.merge(ukp_1,ukp_2,on=["gene_name","uniprot","oid"])
#ukpp["UKBB_PPP_exposure"]=ukpp["gene_name"]+":"+ukpp["uniprot"]+":"+ukpp["panel"]

mr=pd.read_excel("olink_somscan_decode_biobank_protein.xlsx",sheet_name="Sheet3")
mr_decode=mr[~mr["DECODE_exposure"].isna()][["Protein","DECODE_exposure"]].drop_duplicates()
mr_ukbppp=mr[~mr["UKBB_PPP_exposure"].isna()][["Protein","UKBB_PPP_exposure"]].drop_duplicates()
mr_ukbppp["uniprot"]=mr_ukbppp["UKBB_PPP_exposure"].str.split(":",expand=True)[1]
#mr_ukbppp["UKBB_PPP_exposure_1"]=mr_ukbppp["UKBB_PPP_exposure"]


mr_ukbppp_hgnc=pd.merge(mr_ukbppp,hgnc,left_on=["uniprot","Protein"],right_on=["UniProt_ID","Approved symbol"],how="left").drop_duplicates()
mr_ukbppp_not_approved=mr_ukbppp_hgnc[mr_ukbppp_hgnc["UniProt_ID"].isna()].iloc[:,0:3]
mr_ukbppp_approved_symbols=mr_ukbppp_hgnc[mr_ukbppp_hgnc["UniProt_ID"].notna()]

mr_ukbppp_hgnc_previous=pd.merge(mr_ukbppp_not_approved,hgn_previousysmbols,left_on=["uniprot","Protein"],right_on=["UniProt_ID","Previous symbols"],how="left").drop_duplicates()
mr_ukbppp_not_previous=mr_ukbppp_hgnc_previous[mr_ukbppp_hgnc_previous["UniProt_ID"].isna()].iloc[:,0:3]
mr_ukbppp_previous_symbols=mr_ukbppp_hgnc_previous[mr_ukbppp_hgnc_previous["UniProt_ID"].notna()]


mr_ukbppp_hgnc_alias=pd.merge(mr_ukbppp_not_previous,hgnc_alias,left_on=["uniprot","Protein"],right_on=["UniProt_ID",'Alias symbols'],how="left").drop_duplicates()
mr_ukbppp_not_alias=mr_ukbppp_hgnc_alias[mr_ukbppp_hgnc_alias["UniProt_ID"].isna()].iloc[:,0:3]
mr_ukbppp_alias_symbols=mr_ukbppp_hgnc_alias[mr_ukbppp_hgnc_alias["UniProt_ID"].notna()]


##Biomart
mr_ukbppp_biomart_approved=pd.merge(mr_ukbppp_not_alias,biomart_1,left_on=["uniprot","Protein"],right_on=["UniProt_ID",'Approved symbol'],how="left").drop_duplicates()
mr_ukbppp_not_biomart_approved=mr_ukbppp_biomart_approved[mr_ukbppp_biomart_approved["Approved symbol"].isna()].iloc[:,0:3]
mr_ukbppp_biomart_approved=mr_ukbppp_biomart_approved[mr_ukbppp_biomart_approved["Approved symbol"].notna()]
mr_ukbppp_biomart_approved=mr_ukbppp_biomart_approved[mr_ukbppp_biomart_approved["Ensembl_gene_ID"].isin(['ENSG00000079393', 'ENSG00000101307', 'ENSG00000164172','ENSG00000087460',
                                                                               'ENSG00000152061', 'ENSG00000110680','ENSG00000180900'])]

mr_ukbppp_synonym=pd.merge(mr_ukbppp_not_biomart_approved,biomart_2,left_on=["uniprot","Protein"],right_on=["UniProt_ID",'Gene Synonym'],how="left").drop_duplicates()
mr_ukbppp_not_biomart_synonym=mr_ukbppp_synonym[mr_ukbppp_synonym["Approved symbol"].isna()].iloc[:,0:3]
mr_ukbppp_biomart_synonym=mr_ukbppp_synonym[mr_ukbppp_synonym["Approved symbol"].notna()]
mr_ukbppp_not_biomart_synonym[['Ensembl_gene_ID', 'HGNC','Approved symbol', 'NCBI_Gene_ID', 'UniProt_ID', 'OMIM_Gene_ID']]=0


mr_ukpp_withgeneinfo_df=pd.concat([mr_ukbppp_approved_symbols,mr_ukbppp_previous_symbols,mr_ukbppp_alias_symbols,mr_ukbppp_biomart_approved,mr_ukbppp_not_biomart_synonym])
mr_ukpp_withgeneinfo_df=mr_ukpp_withgeneinfo_df.drop(['Previous symbols','Alias symbols','Status'],axis=1)

mr_ukpp_withgeneinfo_df.to_csv("Complete_UKBB-PPP_protein_withcomplete_Geneinformation.csv",index=None)
##aAfter that corrcted uplicate manually

mr_ukpp_withgeneinfo_df=pd.read_csv("Complete_UKBB_PPP_protein_withcomplete_Geneinformation.csv")
mr_ukpp_withgeneinfo_df[['HGNC','OMIM_Gene_ID','NCBI_Gene_ID']]=mr_ukpp_withgeneinfo_df[['HGNC','OMIM_Gene_ID','NCBI_Gene_ID']].astype("Int64")
decode_ukbb_ppp=pd.merge(decode_full,mr_ukpp_withgeneinfo_df,on=['Protein','uniprot', 'HGNC','UniProt_ID', 'OMIM_Gene_ID', 'Ensembl_gene_ID', 'NCBI_Gene_ID',"Approved symbol","Approved name"],how="outer")


decode_full=decode_full.rename(columns={'target_name':'Dcode_target_name','target_full_name':'Dcode_target_full_name',
                                        'platform':'Decode_platform','seqid':'Decode_seqid'})
decode_full=decode_full.drop(['seqid.1', 'seqid2','Status', 'UniProt_ID','Approved name'],axis=1)
decode_full.to_csv("Complete_Deecode_protein_withcomplete_Geneinformation.csv",index=None)


mr_ukpp_withgeneinfo_df2_ukpp=pd.merge(mr_ukpp_withgeneinfo_df2,ukpp,on=["oid","uniprot"],how="left").drop_duplicates()

mr_ukpp_withgeneinfo_df2_ukpp=mr_ukpp_withgeneinfo_df2_ukpp.drop('Approved name',axis=1)

mr_ukpp_withgeneinfo_df2=mr_ukpp_withgeneinfo_df2_ukpp.rename(columns={'oid':'UKBB-PPP_oid','target_full_name':'UKBB-PPP_target_full_name',"gene_name":"UKBB-PPP_target_name",
                                                                  'platform':'UKBB-PPP_platform','panel':'UKBB-PPP_panel'})

mr_ukpp_withgeneinfo_df2_ukpp.to_csv("Complete_UKBB-PPP_protein_withcomplete_Geneinformation.csv",index=None)




### mr_ukpp_withgeneinfo_df2_ukpp
mr_ukpp_withgeneinfo_df2_ukpp=pd.read_csv("Complete_UKBB-PPP_protein_withcomplete_Geneinformation.csv")
mr_ukpp_withgeneinfo_df2_ukpp=mr_ukpp_withgeneinfo_df2_ukpp.rename(columns={'HGNC':'UKBB_PPP_HGNC','Approved symbol':'UKBB_PPP_Approved symbol',
                                'OMIM_Gene_ID':'UKBB_PPP_OMIM_Gene_ID','Ensembl_gene_ID':'UKBB_PPP_Ensembl_gene_ID','NCBI_Gene_ID':'UKBB_PPP_NCBI_Gene_ID'})

### decode_full
decode_full=pd.read_csv("Complete_Decodee_protein_withcomplete_Geneinformation.csv")
decode_full=decode_full.rename(columns={'HGNC':'DCODE_HGNC','Approved symbol':'DCODE_Approved symbol',
                                'OMIM_Gene_ID':'DCODE_OMIM_Gene_ID','Ensembl_gene_ID':'DCODE_Ensembl_gene_ID','NCBI_Gene_ID':'DCODE_NCBI_Gene_ID'})
decode_full=decode_full.rename(columns={'seqid':'DCODE_seqid','target_name':'DCODE_target_name','target_full_name':'DCODE_target_full_name','gene_name':'DCODE_gene_name','platform':'DCODE_platform'})
decode_full=decode_full.drop(["seqid2",'seqid.1'],axis=1)
decode_full=decode_full.drop("UniProt_ID",axis=1)
decode_full=decode_full.rename(columns={'Approved name':'DCODE_Approved name'})



ubb_decode=pd.merge(mr_ukpp_withgeneinfo_df2_ukpp,decode_full,on=['Protein','uniprot'],how="outer")
ubb_decode=ubb_decode.fillna(0)



ubb_decode['HGNC']=np.where(ubb_decode['UKBB_PPP_HGNC']==0,ubb_decode['DCODE_HGNC'],ubb_decode['UKBB_PPP_HGNC'])
ubb_decode['Approved symbol']=np.where(ubb_decode['UKBB_PPP_Approved symbol']==0,ubb_decode['DCODE_Approved symbol'],ubb_decode['UKBB_PPP_Approved symbol'])
ubb_decode['OMIM_Gene_ID']=np.where(ubb_decode['UKBB_PPP_OMIM_Gene_ID']==0,ubb_decode['DCODE_OMIM_Gene_ID'],ubb_decode['UKBB_PPP_OMIM_Gene_ID'])
ubb_decode['Ensembl_gene_ID']=np.where(ubb_decode['UKBB_PPP_Ensembl_gene_ID']==0,ubb_decode['DCODE_Ensembl_gene_ID'],ubb_decode['UKBB_PPP_Ensembl_gene_ID'])
ubb_decode['NCBI_Gene_ID']=np.where(ubb_decode['UKBB_PPP_NCBI_Gene_ID']==0,ubb_decode['DCODE_NCBI_Gene_ID'],ubb_decode['UKBB_PPP_NCBI_Gene_ID'])
ubb_decode['target_full_name']=np.where(ubb_decode['UKBB_PPP_target_full_name']==0,ubb_decode['DCODE_target_full_name'],ubb_decode['UKBB_PPP_target_full_name'])


ubb_decode.to_csv("Complete_UKBB-PPP_Deecode_protein_withcomplete_Geneinformation.csv",index=None)
decode_full.to_csv("Complete_Deecode_protein_withcomplete_Geneinformation.csv",index=None)
mr_ukpp_withgeneinfo_df2_ukpp.to_csv("Complete_UKBB-PPP_protein_withcomplete_Geneinformation.csv",index=None)

ubb_decode_seelected_columns=ubb_decode[['Protein','uniprot','UKBB_PPP_exposure','DECODE_exposure','UKBB_PPP_platform','DCODE_platform',
            'UKBB_PPP_panel','DCODE_seqid','UKBB_PPP_oid','HGNC', 'Approved symbol', 'OMIM_Gene_ID', 'Ensembl_gene_ID','NCBI_Gene_ID', 'target_full_name']]
            
ubb_decode_seelected_columns.to_csv("Complete_UKBB-PPP_Deecode_protein_withcomplete_Geneinformation_SeelecteedColumns.csv",index=None)


########

biomart=biomart[['Ensembl_gene_ID', 'CHR', 'Gene_start', 'Gene_end']].drop_duplicates()

nid=pd.read_csv("teest.csv")
nid.update(nid.select_dtypes(include=['object']).fillna("NA"))
nid['ep_action_mode']=nid['ep_action_mode'].str[:3200]
nid['compound_name']=nid['compound_name'].str[:3200]
nid['compound_id']=nid['compound_id'].str[:3200]

nid.to_csv("Complete_UKBB-PPP_Deecode_protein_withcomplete_Geneinformation_SeelecteedColumns_Annotation_Drug.csv")

test=pd.merge(nid,biomart,on=["Ensembl_gene_ID"],how="left")
test.to_csv("Complete_UKBB-PPP_Deecode_protein_withcomplete_Geneinformation_SeelecteedColumns_Annotation_Drug.csv",index=None)


