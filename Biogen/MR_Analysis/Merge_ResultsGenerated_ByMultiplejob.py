################## Meerge results of cis_exposure_ MR analysis results
import os,glob
import pandas as pd

unique_files=glob.glob("Part0/cis_exposure_*")
unique_files=[x.split("/")[1] for x in unique_files ]

print("\n\n",len(unique_files))
for i in unique_files:
    print(i)

for file in unique_files:
    Files=glob.glob(f"Part*/{file}")
    master_df=pd.DataFrame()
    for part_file in Files:
        part_file_df=pd.read_csv(part_file)
        master_df=pd.concat([master_df,part_file_df])
    master_df.to_csv(f"{file}",index=None)





## trans_exposure_TwoSampleMR_
unique_files=glob.glob("Part0/trans_exposure_*")
unique_files=[x.split("/")[1] for x in unique_files ]
unique_files=[x for x in unique_files if "MHC" not in x]

print("\n\n",len(unique_files))
for i in unique_files:
    print(i)
    
for file in unique_files:
    Files=glob.glob(f"Part*/{file}")
    master_df=pd.DataFrame()
    for part_file in Files:
        part_file_df=pd.read_csv(part_file)
        master_df=pd.concat([master_df,part_file_df])
    master_df.to_csv(f"{file}",index=None)



##trans_exposure_noMHC_
unique_files=glob.glob("Part0/trans_exposure_noMHC_*")
unique_files=[x.split("/")[1] for x in unique_files ]

print("\n\n",len(unique_files))
for i in unique_files:
    print(i)
    
for file in unique_files:
    Files=glob.glob(f"Part*/{file}")
    master_df=pd.DataFrame()
    for part_file in Files:
        part_file_df=pd.read_csv(part_file)
        master_df=pd.concat([master_df,part_file_df])
    master_df.to_csv(f"{file}",index=None)
    




##trans_exposure_NoMHC_Unique_
unique_files=glob.glob("Part0/trans_exposure_NoMHC_Unique_*")
unique_files=[x.split("/")[1] for x in unique_files ]

print("\n\n",len(unique_files))
for i in unique_files:
    print(i)

for file in unique_files:
    Files=glob.glob(f"Part*/{file}")
    master_df=pd.DataFrame()
    for part_file in Files:
        part_file_df=pd.read_csv(part_file)
        master_df=pd.concat([master_df,part_file_df])
    master_df.to_csv(f"{file}",index=None)



os.system("mkdir -p Compiled_MR_Results/cis_exposure")
os.system("mkdir -p Compiled_MR_Results/trans_exposure")
os.system("mkdir -p Compiled_MR_Results/trans_exposure_noMHC")
os.system("mkdir -p Compiled_MR_Results/trans_exposure_NoMHC_Unique")

os.system("mv trans_exposure_NoMHC_Unique_* Compiled_MR_Results/trans_exposure_NoMHC_Unique")
os.system("mv trans_exposure_noMHC_* Compiled_MR_Results/trans_exposure_noMHC")
os.system("mv trans_exposure_* Compiled_MR_Results/trans_exposure")
os.system("mv cis_exposure_* Compiled_MR_Results/cis_exposure")

