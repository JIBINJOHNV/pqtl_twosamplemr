import pandas as pd

vcf=pd.read_csv("Decode_Proteein.txt",header=None)
vcf[0]=vcf[0].str.split("/",expand=True)[6]
vcf["Symbol"]=vcf[0].str.split("_",expand=True)[2]
vcf["Sequence_ID"]=vcf[0].str.split("_",expand=True)[0]+"-"+vcf[0].str.split("_",expand=True)[1]
soma=pd.read_csv("SomaScanMenu.csv",sep=",")
vcf_soma=pd.merge(vcf,soma,left_on=["Sequence_ID","Symbol"],right_on=["Sequence_ID","EntrezeeGene_Nname"],how="left")
vcf_soma_missing=vcf_soma[vcf_soma["Target_Name"].isna()].iloc[:,0:3]
vcf_soma_notmissing=vcf_soma[~vcf_soma["Target_Name"].isna()]



mart=pd.read_csv("mart_export.txt")
mart2=mart[['Gene name',"UniProtKB/Swiss-Prot ID","Chromosome/scaffold name",'Transcript start (bp)', 'Transcript end (bp)']].drop_duplicates()
mart3=mart2.groupby(["Gene name","UniProtKB/Swiss-Prot ID","Chromosome/scaffold name"]).agg({"Transcript start (bp)": "min","Transcript end (bp)": "max"}).reset_index()
mart3_haploid=mart3[mart3["Chromosome/scaffold name"].str.contains("H",na=False)]
mart3_nothaploid=mart3[~mart3["Chromosome/scaffold name"].str.contains("H",na=False)]
mart3_haploid_not_inmart3_nothaploid=mart3_haploid[~mart3_haploid["Gene name"].isin(mart3_nothaploid["Gene name"].unique())]

mergeeed=pd.merge(vcf_soma_notmissing,mart3_nothaploid,left_on=["UniprotID","Symbol"],right_on=["UniProtKB/Swiss-Prot ID","Gene name"],how="left")
merge_missing=mergeeed[mergeeed["Chromosome/scaffold name"].isna()].iloc[:,0:4]
merge_nomissing=mergeeed[~mergeeed["Chromosome/scaffold name"].isna()]


missing=pd.concat([merge_missing,vcf_soma_missing]).reset_index().drop("index",axis=1)

##Some case uniprot id ot matching , so only merge with Gene name
mart=pd.read_csv("mart_export.txt")
mart2=mart[['Gene name',"Chromosome/scaffold name",'Transcript start (bp)', 'Transcript end (bp)']].drop_duplicates()
mart3=mart2.groupby(["Gene name","Chromosome/scaffold name"]).agg({"Transcript start (bp)": "min","Transcript end (bp)": "max"}).reset_index()
mart3_haploid=mart3[mart3["Chromosome/scaffold name"].str.contains("H",na=False)]
mart3_nothaploid=mart3[~mart3["Chromosome/scaffold name"].str.contains("H",na=False)]
mart3_haploid_not_inmart3_nothaploid=mart3_haploid[~mart3_haploid["Gene name"].isin(mart3_nothaploid["Gene name"].unique())]

mergeeed2=pd.merge(missing,mart3_nothaploid,left_on="Symbol",right_on="Gene name",how="left")
merge_missing2=mergeeed2[mergeeed2["Chromosome/scaffold name"].isna()].iloc[:,0:4]
merge_nomissing2=mergeeed2[~mergeeed2["Chromosome/scaffold name"].isna()]


##Some case uniprot id and name not match so Synoms used for matching
mart=pd.read_csv("mart_export.txt")
mart2=mart[['Gene Synonym',"Chromosome/scaffold name",'Transcript start (bp)', 'Transcript end (bp)']].drop_duplicates()
mart3=mart2.groupby(["Gene Synonym","Chromosome/scaffold name"]).agg({"Transcript start (bp)": "min","Transcript end (bp)": "max"}).reset_index()
mart3_haploid=mart3[mart3["Chromosome/scaffold name"].str.contains("H",na=False)]
mart3_nothaploid=mart3[~mart3["Chromosome/scaffold name"].str.contains("H",na=False)]
mart3_haploid_not_inmart3_nothaploid=mart3_haploid[~mart3_haploid["Gene Synonym"].isin(mart3_nothaploid["Gene Synonym"].unique())]

mergeeed3=pd.merge(merge_missing2,mart3_nothaploid,left_on="Symbol",right_on="Gene Synonym",how="left")
merge_missing3=mergeeed3[mergeeed3["Chromosome/scaffold name"].isna()].iloc[:,0:4]
merge_nomissing3=mergeeed3[~mergeeed3["Chromosome/scaffold name"].isna()]
merge_nomissing3=pd.merge(merge_nomissing3,mart[["Gene name", "Gene Synonym"]].drop_duplicates(),on="Gene Synonym")

uniq=merge_nomissing3[~merge_nomissing3[0].duplicated(keep=False)]
duplicated=merge_nomissing3[merge_nomissing3[0].duplicated(keep=False)]
duplicated["Chromosome/scaffold name"]=duplicated["Chromosome/scaffold name"].astype("str")

duplicated_uniq=pd.DataFrame()
duplicated_uniq=pd.concat([duplicated_uniq,duplicated[(duplicated["Chromosome/scaffold name"]=="7")&(duplicated["Gene name"]=="PALS2" )]])
duplicated_uniq=pd.concat([duplicated_uniq,duplicated[(duplicated["Chromosome/scaffold name"]=="12")&(duplicated["Gene name"]=="KMT2B" )]])
duplicated_uniq=pd.concat([duplicated_uniq,duplicated[(duplicated["Chromosome/scaffold name"]=="18")&(duplicated["Gene name"]=="ADCYAP1" )]])
duplicated_uniq=pd.concat([duplicated_uniq,duplicated[(duplicated["Chromosome/scaffold name"]=="11")&(duplicated["Gene name"]=="CBLIF" )]])
duplicated_uniq=pd.concat([duplicated_uniq,duplicated[(duplicated["Chromosome/scaffold name"]=="1")&(duplicated["Gene name"]=="SARS1" )]])
duplicated_uniq=pd.concat([duplicated_uniq,duplicated[(duplicated["Chromosome/scaffold name"]=="10")&(duplicated["Gene name"]=="CDK1" )]])
duplicated_uniq=pd.concat([duplicated_uniq,duplicated[(duplicated["Chromosome/scaffold name"]=="8")&(duplicated["Gene name"]=="CCN3" )]])

merge_nomissing3=pd.concat([uniq,duplicated_uniq])

Result_Notmissing=pd.concat([merge_nomissing,merge_nomissing2,merge_nomissing3])
Result_Notmissing2=Result_Notmissing[[0,'Symbol','Sequence_ID','Target_Name','EntrezeeGene_Nname','Gene name','UniprotID','UniProtKB/Swiss-Prot ID',
                      'Chromosome/scaffold name','Transcript start (bp)','Transcript end (bp)']]

Result_Notmissing3=Result_Notmissing2[[0,'Symbol','Sequence_ID','Gene name','Chromosome/scaffold name','Transcript start (bp)','Transcript end (bp)']]
Result_Notmissing3.columns=["Filename","SomaScan_Symbol",'SomaScan_Sequence_ID','Gene name','Biomart_Chromosome','Biomart_TranscriptStart','Biomart_TranscriptEnd']
Result_Notmissing3[["Biomart_TranscriptStart",  "Biomart_TranscriptEnd"]]=Result_Notmissing3[["Biomart_TranscriptStart",  "Biomart_TranscriptEnd"]].astype("int")


Result_Notmissing3.to_csv("Decode_Proteins_TrnnscriptSTart_End_Hg38.csv",index=None)
merge_missing3.to_csv("Decode_Proteins_DontHave_TrnnscriptSTart_End_Hg38.csv",index=None)
