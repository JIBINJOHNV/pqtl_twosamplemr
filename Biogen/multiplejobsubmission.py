import pandas as pd
import numpy as np
import os


# Load your original DataFrame
df = pd.read_csv("/edgehpc/dept/human_genetics/users/jjohn1/Cis_Trans/Exposure/UKB_50KClumped/Biogen/Reanalysis/All_Biogen_Exposure_TSS_chr_index_3k.csv")


# Split the DataFrame into almost equal-sized parts
split_dfs = np.array_split(df, 50)

basedir=os.getcwd()+"/"

# Save each part to separate files
for i, split_df in enumerate(split_dfs):
    os.system(f"mkdir -p Part{i}")
    Dir=f"{basedir}Part{i}/"
    split_df.to_csv(f"{Dir}All_Biogen_Exposure_TSS_chr_index_3k_Part{i}.csv", index=False)
    os.system(f"cp sbatch_command.sh {Dir}")
    os.system(f"sed s'|All_Biogen_Exposure_TSS_chr_index_3k_Part1.csv|{Dir}All_Biogen_Exposure_TSS_chr_index_3k_Part{i}.csv|' pQTL_Expsure_selection_clumping.R > {basedir}Part{i}/pQTL_Expsure_selection_clumping.R")
    os.chdir(Dir)
    os.system(f"sbatch sbatch_command.sh ")
    os.chdir(basedir)









################## Meerge results from above job
import os,glob
import pandas as pd

Files=glob.glob("Part*/All_significannt_CisTransexposure_BeforeQC_LDclumping.csv")
Folders=sorted([x.replace("/All_significannt_CisTransexposure_BeforeQC_LDclumping.csv","") for x in Files])
Folders=[x for x in Folders  if x!='Part1']
basedir=os.getcwd()+"/"


for foder in Folders:
    Dir=f"{basedir}{foder}/"
    os.system(f"cp Decode_twosampleMR_Running.R mr_analysis_sbatch_command.sh {Dir}")
    os.chdir(Dir)
    os.system(f"sbatch mr_analysis_sbatch_command.sh ")
    os.chdir(basedir)
