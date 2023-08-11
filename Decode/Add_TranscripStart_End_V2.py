from pybiomart import Server
import pandas as pd
import os

os.system("wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00978-w/MediaObjects/41588_2021_978_MOESM4_ESM.xlsx") #From Decode publication

server = Server(host='http://www.ensembl.org')

dataset = (server.marts['ENSEMBL_MART_ENSEMBL']
                 .datasets['hsapiens_gene_ensembl'])

atributeees=dataset.list_attributes()

res1=dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','chromosome_name',"transcription_start_site","transcript_end","transcript_start","entrezgene_id","hgnc_id"])
res2=dataset.query(attributes=['ensembl_gene_id', 'external_gene_name',"uniprotswissprot","uniprotsptrembl"])
mart=pd.merge(res1,res2,on=['Gene stable ID', 'Gene name'],how="outer")
mart2=mart[['Gene name',"UniProtKB/Swiss-Prot ID","Chromosome/scaffold name",'Transcript start (bp)', 'Transcript end (bp)','Gene stable ID',"NCBI gene (formerly Entrezgene) ID","HGNC ID"]].drop_duplicates()
mart2.rename(columns={'Gene name':"Biomart_GeneName", 'UniProtKB/Swiss-Prot ID':'Biomart_UniProtKB/Swiss-ProtID','Gene stable ID':'Biomart_GeneStableID','NCBI gene (formerly Entrezgene) ID':'Biomart_NCBI_EntrezgeneID', 'HGNC ID':'BiomartHGNCID'},inplace=True)
mart2=mart2.drop_duplicates()
mart2["Chromosome/scaffold name"]=mart2["Chromosome/scaffold name"].astype("str")
mart2_Nochr=mart2[mart2["Chromosome/scaffold name"].str.contains("H|G|K")]
mart2=mart2[~mart2["Chromosome/scaffold name"].str.contains("H|G|K")]

decode_pub=pd.read_excel("41588_2021_978_MOESM4_ESM.xlsx",sheet_name='ST01')
decode_pub=pd.read_excel("41588_2021_978_MOESM4_ESM.xlsx",sheet_name='ST01',skiprows=[0,1])
decode_pub2=decode_pub[['SeqId', 'Protein (short name)', 'Protein (full name)', 'Gene','UniProt','Ensembl.Gene.ID', 'Entrez.Gene.ID', 'HGNC.ID']].drop_duplicates()
decode_pub2.rename(columns={'SeqId':"Sequence_ID", 'Protein (short name)':'DecodePub_Protein(short name)', 'Protein (full name)':'DecodePub_Protein(full name)', 'Gene':'DecodePub_Gene',
                'UniProt':'DecodePub_UniProt', 'Ensembl.Gene.ID':'DecodePub_Ensembl.Gene.ID', 'Entrez.Gene.ID':'DecodePub_Entrez.Gene.ID', 'HGNC.ID':'DecodePub_HGNC.ID' },inplace=True)


vcf=pd.read_csv("Decode_Proteein.txt",header=None)
vcf[0]=vcf[0].str.split("/",expand=True)[6]
vcf["Symbol"]=vcf[0].str.split("_",expand=True)[2]
#vcf["Symbol2"]=vcf[0].str.split("_",expand=True)[3].str.replace(".txt.gz","")
vcf["Sequence_ID"]=vcf[0].str.split("_",expand=True)[0]+"-"+vcf[0].str.split("_",expand=True)[1]
vcf["Sequence_ID"]=vcf["Sequence_ID"].str.replace("-","_")

decode_pub2_vcf=pd.merge(vcf,decode_pub2,on="Sequence_ID")





decode_pub2_vcf_mart2_1=pd.merge(decode_pub2_vcf,mart2,left_on=['DecodePub_UniProt','DecodePub_Ensembl.Gene.ID','DecodePub_Entrez.Gene.ID','DecodePub_HGNC.ID'], 
                                                      right_on=['Biomart_UniProtKB/Swiss-ProtID','Biomart_GeneStableID','Biomart_NCBI_EntrezgeneID', 'BiomartHGNCID'],how="left")
decode_pub2_vcf_mart2_Notmssing_1=decode_pub2_vcf_mart2_1[~decode_pub2_vcf_mart2_1["Transcript start (bp)"].isna()]
decode_pub2_vcf_mart2_mssing_1=decode_pub2_vcf_mart2_1[decode_pub2_vcf_mart2_1["Transcript start (bp)"].isna()].iloc[:,:-8]

decode_pub2_vcf_mart2_2=pd.merge(decode_pub2_vcf_mart2_mssing_1,mart2,left_on=['DecodePub_UniProt','DecodePub_Entrez.Gene.ID','DecodePub_HGNC.ID'], 
                                                      right_on=['Biomart_UniProtKB/Swiss-ProtID','Biomart_NCBI_EntrezgeneID', 'BiomartHGNCID'],how="left")
decode_pub2_vcf_mart2_Notmssing_2=decode_pub2_vcf_mart2_2[~decode_pub2_vcf_mart2_2["Transcript start (bp)"].isna()]
decode_pub2_vcf_mart2_mssing_2=decode_pub2_vcf_mart2_2[decode_pub2_vcf_mart2_2["Transcript start (bp)"].isna()].iloc[:,:-8]
                                                      
decode_pub2_vcf_mart2_3=pd.merge(decode_pub2_vcf_mart2_mssing_2,mart2,left_on=['DecodePub_UniProt','DecodePub_HGNC.ID'], 
                                                      right_on=['Biomart_UniProtKB/Swiss-ProtID', 'BiomartHGNCID'],how="left")
decode_pub2_vcf_mart2_Notmssing_3=decode_pub2_vcf_mart2_3[~decode_pub2_vcf_mart2_3["Transcript start (bp)"].isna()]
decode_pub2_vcf_mart2_mssing_3=decode_pub2_vcf_mart2_3[decode_pub2_vcf_mart2_3["Transcript start (bp)"].isna()].iloc[:,:-8]
                                                 
decode_pub2_vcf_mart2_4=pd.merge(decode_pub2_vcf_mart2_mssing_3,mart2,left_on=['DecodePub_UniProt'], right_on=['Biomart_UniProtKB/Swiss-ProtID'],how="left")
decode_pub2_vcf_mart2_Notmssing_4=decode_pub2_vcf_mart2_4[~decode_pub2_vcf_mart2_4["Transcript start (bp)"].isna()]
decode_pub2_vcf_mart2_mssing_4=decode_pub2_vcf_mart2_4[decode_pub2_vcf_mart2_4["Transcript start (bp)"].isna()].iloc[:,:-8]
                                      
decode_pub2_vcf_mart2_5=pd.merge(decode_pub2_vcf_mart2_mssing_4,mart2,left_on=['DecodePub_HGNC.ID'], right_on=['BiomartHGNCID'],how="left")
decode_pub2_vcf_mart2_Notmssing_5=decode_pub2_vcf_mart2_5[~decode_pub2_vcf_mart2_5["Transcript start (bp)"].isna()]
decode_pub2_vcf_mart2_mssing_5=decode_pub2_vcf_mart2_5[decode_pub2_vcf_mart2_5["Transcript start (bp)"].isna()].iloc[:,:-8]
                                                                                                    
decode_pub2_vcf_mart2_6=pd.merge(decode_pub2_vcf_mart2_mssing_5,mart2,left_on=['DecodePub_Entrez.Gene.ID'], right_on=['Biomart_NCBI_EntrezgeneID'],how="left")
decode_pub2_vcf_mart2_Notmssing_6=decode_pub2_vcf_mart2_6[~decode_pub2_vcf_mart2_6["Transcript start (bp)"].isna()]
decode_pub2_vcf_mart2_mssing_6=decode_pub2_vcf_mart2_6[decode_pub2_vcf_mart2_6["Transcript start (bp)"].isna()].iloc[:,:-8]
                                                                        
decode_pub2_vcf_mart2_7=pd.merge(decode_pub2_vcf_mart2_mssing_6,mart2,left_on=['DecodePub_Ensembl.Gene.ID'], right_on=['Biomart_GeneStableID'],how="left")
decode_pub2_vcf_mart2_Notmssing_7=decode_pub2_vcf_mart2_7[~decode_pub2_vcf_mart2_7["Transcript start (bp)"].isna()]
decode_pub2_vcf_mart2_mssing_7=decode_pub2_vcf_mart2_7[decode_pub2_vcf_mart2_7["Transcript start (bp)"].isna()].iloc[:,:-8]


##No chr
decode_pub2_vcf_mart2_8=pd.merge(decode_pub2_vcf_mart2_mssing_7,mart2_Nochr,left_on=['DecodePub_UniProt'], right_on=['Biomart_UniProtKB/Swiss-ProtID'],how="left")
decode_pub2_vcf_mart2_Notmssing_8=decode_pub2_vcf_mart2_8[~decode_pub2_vcf_mart2_8["Transcript start (bp)"].isna()]
decode_pub2_vcf_mart2_mssing_8=decode_pub2_vcf_mart2_8[decode_pub2_vcf_mart2_8["Transcript start (bp)"].isna()].iloc[:,:-8]
decode_pub2_vcf_mart2_Notmssing_8=decode_pub2_vcf_mart2_Notmssing_8.drop_duplicates(subset=[0,"Sequence_ID"])

##For this genee no chromosome position was available in Ensembl or NCBI, but found im OMIM, SO USED THAT POSITION https://www.omim.org/entry/617777
BAGE3_Dict={'Biomart_GeneName':"BAGE3", 'Biomart_UniProtKB/Swiss-ProtID':"Q86Y29",'Chromosome/scaffold name':"21", 'Transcript start (bp)':7000001,
'Transcript end (bp)':10900000, 'Biomart_GeneStableID':"NA",'Biomart_NCBI_EntrezgeneID':85318, 'BiomartHGNCID':"HGNC:15728"}
BAGE3_series=pd.Series(BAGE3_Dict)
BAGE3_df=pd.concat([pd.DataFrame(BAGE3_series),decode_pub2_vcf_mart2_mssing_8.T],axis=0)
BAGE3_df[0]=BAGE3_df[0]=np.where(BAGE3_df[0].isna(),BAGE3_df[50],BAGE3_df[0])
BAGE3_df=BAGE3_df[[0]].T



decode_pub2_vcf_mart_merged=pd.concat([decode_pub2_vcf_mart2_Notmssing_1,decode_pub2_vcf_mart2_Notmssing_2,decode_pub2_vcf_mart2_Notmssing_3,
                                       decode_pub2_vcf_mart2_Notmssing_4,decode_pub2_vcf_mart2_Notmssing_5,decode_pub2_vcf_mart2_Notmssing_6,
                                       decode_pub2_vcf_mart2_Notmssing_7,decode_pub2_vcf_mart2_Notmssing_8,BAGE3_df])



decode_pub2_vcf_mart_merged_2=decode_pub2_vcf_mart_merged[[0,'Symbol','Sequence_ID','Chromosome/scaffold name','Transcript start (bp)','Transcript end (bp)']].drop_duplicates()
mart3=decode_pub2_vcf_mart_merged_2.groupby([0,"Symbol","Sequence_ID","Chromosome/scaffold name"]).agg({"Transcript start (bp)": "min","Transcript end (bp)": "max"}).reset_index()
decode_pub2_vcf_mart=pd.merge(mart3,decode_pub2_vcf,on=[0,"Symbol", "Sequence_ID"])



decode_pub2_vcf_mart2=decode_pub2_vcf_mart[[0, 'DecodePub_Gene','Sequence_ID','Symbol','Chromosome/scaffold name','Transcript start (bp)','Transcript end (bp)']]
decode_pub2_vcf_mart2.columns=["Filename","SomaScan_Symbol",'SomaScan_Sequence_ID','Gene name','Biomart_Chromosome','Biomart_TranscriptStart','Biomart_TranscriptEnd']
decode_pub2_vcf_mart2[["Biomart_TranscriptStart",  "Biomart_TranscriptEnd"]]=decode_pub2_vcf_mart2[["Biomart_TranscriptStart",  "Biomart_TranscriptEnd"]].astype("int")



decode_pub2_vcf_mart2.to_csv("Decode_Proteins_TrnnscriptSTart_End_Hg38.csv",index=None)
decode_pub2_vcf_mart.to_csv("Decode_Proteins_TrnnscriptSTart_End_Hg38_WithCompleteInformation.csv",index=None)

