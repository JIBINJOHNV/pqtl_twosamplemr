
import pandas as pd

gene_info=pd.read_csv("Decode_Biogen_CisExposure_TransExposureNoMHC_suppelementary_table_withGeeneInfo.csv")
gene_info=gene_info[['UniProt_ID', 'Approved symbol', 'HGNC', 'Ensembl_gene_ID', 'NCBI_Gene_ID']].drop_duplicates()
gene_info=gene_info.rename(columns={"UniProt_ID":"uniprot_ID"})
gene_info=gene_info[gene_info['Ensembl_gene_ID']!="ENSG00000278145"].drop_duplicates()

gene_info_list=[{"uniprot_ID":"Q8WXW3-4","HGNC":23352,"Ensembl_gene_ID":"ENSG00000083535","NCBI_Gene_ID":10464,"Approved symbol":'PIBF1'},
 {"uniprot_ID":"P0DN86","HGNC":1886,"Ensembl_gene_ID":"ENSG00000104827","NCBI_Gene_ID":1082,"Approved symbol":'CGB3'},
 {"uniprot_ID":"Q9NV35","HGNC":23063,"Ensembl_gene_ID":"ENSG00000136159","NCBI_Gene_ID":55270,"Approved symbol":'NUDT15' },
 {"uniprot_ID":"Q8N6C8","HGNC":6604,"Ensembl_gene_ID":"ENSG00000275841","NCBI_Gene_ID":11026,"Approved symbol":'LILRA3' },
 {"uniprot_ID":"O75882-2","HGNC":885,"Ensembl_gene_ID":"ENSG00000088812","NCBI_Gene_ID":8455,"Approved symbol":'ATRN' },
 {"uniprot_ID":"NTproBNP","HGNC":7940,"Ensembl_gene_ID":"ENSG00000120937","NCBI_Gene_ID":4879,"Approved symbol":'NPPB' },
 {"uniprot_ID":"Q29980_Q29983","Approved symbol":'MICA_MICB' },
 {"uniprot_ID":"P21217_Q11128","Approved symbol":'FUT3_FUT5' },
 {"uniprot_ID":"P0DUB6_P0DTE7_P0DTE8","Approved symbol":'AMY1A' },
 {"uniprot_ID":"Q14213_Q8NEV9","Approved symbol":'EBI3_IL27' }]

gene_info_list_df=pd.DataFrame(gene_info_list)
gene_info=pd.concat([gene_info,gene_info_list_df])



biogen_cis=pd.read_csv("Biogen_All_significannt_cis_exposure_AfterQC_LDclumping.csv")
biogen_cis['uniprot']=biogen_cis["exposure"].str.split(":",expand=True)[1]
biogen_cis=biogen_cis.drop("Unnamed: 0",axis=1)
biogen_cis['panel']=biogen_cis["exposure"].str.split(":",expand=True)[2]
biogen_cis['panel']=biogen_cis['panel'].str.split("_",expand=True)[0]

biogen_details_1=pd.read_csv("ukbb-ppp_protein_details.txt",sep=" ",header=None)
biogen_details_1.columns=["gene_name","uniprot","oid","v2","panel"]
biogen_details_1["platform"]="Olink Explore"

biogen_cis_with=pd.merge(biogen_cis,biogen_details_1,left_on=["id","uniprot",'panel'],right_on=["gene_name","uniprot",'panel'],how="left").drop_duplicates()
biogen_cis_with_oid=biogen_cis_with[biogen_cis_with["gene_name"].notna()].drop_duplicates()
biogen_cis_with_no_oid=biogen_cis_with[biogen_cis_with["gene_name"].isna()].drop(['gene_name','oid','v2','platform'],axis=1).drop_duplicates()

biogen_details=pd.read_excel("olink_somscan_decode_biobank_protein.xlsx",sheet_name="Sheet2")
biogen_details=biogen_details[['gene_name',"oid","uniprot",'platform']].drop_duplicates()
custom_olinkanno=[
                    {"uniprot":'P47929',"gene_name":'LGALS7_LGALS7B',"oid":'OID21406'},
                    {"uniprot":'P0DUB6_P0DTE7_P0DTE8',"gene_name":'AMY1A',"oid":'OID30707'},
                    {"uniprot":'Q9H3K6',"gene_name":'BOLA2',"oid":'OID30267'},
                    {"uniprot":'Q96QH8',"gene_name":'SPACA5',"oid":'OID30990'},
                    {"uniprot":'P0DN86',"gene_name":'CGB3',"oid":'OID30671'},
                    {"uniprot":'O15263',"gene_name":'DEFB4A_DEFB4B',"oid":'OID21373'},
                    {"uniprot":'Q14213_Q8NEV9',"gene_name":'EBI3_IL27',"oid":'OID21389'},
                    {"uniprot":'Q29980_Q29983',"gene_name":'MICA_MICB',"oid":'OID20593'},
                    {"uniprot":'P21217_Q11128',"gene_name":'FUT3_FUT5',"oid":'OID21013'},
                    {"uniprot":'P59665',"gene_name":'DEFA1_DEFA1B',"oid":'OID20344'},
                    {"uniprot":'Q8WTQ1',"gene_name":'DEFB104A_DEFB104B',"oid":'OID30849'},
                    {"uniprot":'P12532',"gene_name":'CKMT1A_CKMT1B',"oid":'OID20721'}]

biogen_details_2=pd.DataFrame(custom_olinkanno)
biogen_details_2['platform']="Olink Explore"
biogen_cis_with_oid_2=pd.merge(biogen_cis_with_no_oid,biogen_details_2,on=["uniprot"],how="left")
biogen_cis=pd.concat([biogen_cis_with_oid_2,biogen_cis_with_oid]).drop("v2",axis=1)






decode_cis=pd.read_csv("Decode_All_significannt_cis_exposure_AfterQC_LDclumping.csv")
decode_details=pd.read_excel("olink_somscan_decode_biobank_protein.xlsx",sheet_name="Sheet1")
decode_details=decode_details[["seqid2","uniprot",'platform']].drop_duplicates()
decode_cis=pd.merge(decode_cis,decode_details,left_on="exposure",right_on="seqid2",how="left").drop("seqid2",axis=1)
decode_cis["seqid"]=decode_cis["exposure"].str.split("_",expand=True)[1]+"_"+decode_cis["exposure"].str.split("_",expand=True)[2]



cis_merged=pd.merge(decode_cis[["uniprot","exposure"]].drop_duplicates(),biogen_cis[["uniprot",
                                                "exposure"]].drop_duplicates(),on="uniprot",how="outer")


cis_merged[ (cis_merged["exposure_y"].notna() & cis_merged["exposure_x"].notna()) ]
cis_merged[ (cis_merged["exposure_y"].notna() & cis_merged["exposure_x"].isna()) ]
cis_merged[ (cis_merged["exposure_x"].notna() & cis_merged["exposure_y"].isna()) ]

len(cis_merged[ (cis_merged["exposure_y"].notna() & cis_merged["exposure_x"].notna()) ]['uniprot'].unique())
len(cis_merged[ (cis_merged["exposure_y"].notna() & cis_merged["exposure_x"].isna()) ]['uniprot'].unique())
len(cis_merged[ (cis_merged["exposure_x"].notna() & cis_merged["exposure_y"].isna()) ]['uniprot'].unique())



biogen_cis.columns=[x+"_UKBB-PPP" for x in biogen_cis.columns ]
biogen_cis=biogen_cis.rename(columns={"uniprot_UKBB-PPP":"uniprot_ID","ID_UKBB-PPP":"Variant_ID",'seqnames_UKBB-PPP':'seqnames','start_UKBB-PPP':'start','REF_UKBB-PPP':'REF','ALT_UKBB-PPP':'ALT'})

decode_cis.columns=[x+"_Decode" for x in decode_cis.columns ]
decode_cis=decode_cis.rename(columns={"uniprot_Decode":"uniprot_ID","ID_Decode":"Variant_ID",'seqnames_Decode':'seqnames','start_Decode':'start',"REF_Decode":"REF","ALT_Decode":"ALT"})

cis_df=pd.merge(biogen_cis,decode_cis,on=["Variant_ID","uniprot_ID",'seqnames','start','REF','ALT'],how="outer")
cis_df_annot=pd.merge(cis_df,gene_info,on="uniprot_ID",how="left")
cis_df_annot=cis_df_annot[['uniprot_ID','Approved symbol', 'HGNC', 'Ensembl_gene_ID', 'NCBI_Gene_ID','Variant_ID','seqnames', 'start', 'REF', 'ALT', 'id_UKBB-PPP', 'ES_UKBB-PPP', 'SE_UKBB-PPP', 'P_UKBB-PPP', 'AF_UKBB-PPP', 'SI_UKBB-PPP', 'SS_UKBB-PPP', 'exposure_UKBB-PPP','oid_UKBB-PPP', 'platform_UKBB-PPP', 'id_Decode','ES_Decode', 'SE_Decode', 'P_Decode', 'AF_Decode', 'SI_Decode', 'SS_Decode', 'exposure_Decode', 'platform_Decode', 'seqid_Decode']]
cis_df_annot.to_csv("UKBB-PPP_Decode_All_significannt_cis_exposure_AfterQC_LDclumping_Pqtl_level.csv",index=None)




biogen_cis2=biogen_cis[['uniprot_ID','id_UKBB-PPP', 'exposure_UKBB-PPP','oid_UKBB-PPP', 'platform_UKBB-PPP']].drop_duplicates()
decode_cis2=decode_cis[['uniprot_ID','id_Decode','exposure_Decode', 'platform_Decode', 'seqid_Decode']].drop_duplicates()
cis_df2=pd.merge(biogen_cis2,decode_cis2,on=["uniprot_ID"],how="outer").drop_duplicates()
cis_df_annot2=pd.merge(cis_df2,gene_info,on="uniprot_ID",how="left").drop_duplicates()
cis_df_annot2=cis_df_annot2[['uniprot_ID','Approved symbol', 'HGNC', 'Ensembl_gene_ID','NCBI_Gene_ID','id_UKBB-PPP','exposure_UKBB-PPP','oid_UKBB-PPP','platform_UKBB-PPP','id_Decode', 'exposure_Decode', 'platform_Decode','seqid_Decode']]
cis_df_annot2.to_csv("UKBB-PPP_Decode_All_significannt_cis_exposure_AfterQC_LDclumping_Protein_level.csv",index=None)


