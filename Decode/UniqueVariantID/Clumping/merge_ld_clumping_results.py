
import pandas as pd
import numpy as np
import os,glob

Files=["All_significannt_cis_exposure_AfterQC_LDclumping.csv","All_significannt_trans_exposure_AfterQC_LDclumping.csv",
       "All_significannt_trans_exposure_AfterQC_LDclumping_MHCRemoval.csv","All_significannt_CisTransexposure_BeforeQC_LDclumping.csv"]

os.system(f"mkdir CombinedResultsfromAllBatches")

for file in Files:
    specific_files=glob.glob(f"Part*/{file}")
    master_df=pd.DataFrame()
    for specific_file in specific_files:
        temp_df=pd.read_csv(specific_file)
        master_df=pd.concat([master_df,temp_df])
    
    master_df.to_csv(f"CombinedResultsfromAllBatches/Decode_{file}",index=None)


##Creeatee uniq mhc free trans variant    
trans_nomhc=pd.read_csv(f"CombinedResultsfromAllBatches/Decode_All_significannt_trans_exposure_AfterQC_LDclumping_MHCRemoval.csv")
trans_nomhc.drop("Unnamed: 0",axis=1,inplace=True)
trans_nomhc_count=trans_nomhc.groupby("ID")["seqnames"].count().reset_index().rename(columns={"seqnames":"Count"})

merged=pd.merge(trans_nomhc,trans_nomhc_count,on="ID")
merged_unique=merged[merged["Count"]==1]
merged_unique=merged_unique.drop_duplicates().drop("Count",axis=1)
merged_unique.to_csv(f"CombinedResultsfromAllBatches/Decode_All_significannt_trans_exposure_AfterQC_LDclumping_NoMHC_Unique.csv",index=None)

