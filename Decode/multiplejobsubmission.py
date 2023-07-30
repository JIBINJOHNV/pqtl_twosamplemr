
import pandas as pd
import numpy as np
import os


# Load your original DataFrame
df = pd.read_csv("Decode_Proteins_TrnnscriptSTart_End_Hg38.csv")


# Split the DataFrame into almost equal-sized parts
split_dfs = np.array_split(df, 50)

basedir=os.getcwd()+"/"

# Save each part to separate files
for i, split_df in enumerate(split_dfs):
    os.system(f"mkdir -p Part{i}")
    Dir=f"{basedir}Part{i}/"
    split_df.to_csv(f"{Dir}Decode_Proteins_TrnnscriptSTart_End_Hg38_Part{i}.csv", index=False)
    os.system(f"cp sbatch_command.sh {Dir}")
    os.system(f"sed s'|Decode_Proteins_TrnnscriptSTart_End_Hg38_Part1.csv|{Dir}Decode_Proteins_TrnnscriptSTart_End_Hg38_Part{i}.csv|' pQTL_Expsure_selection_clumping.R > {basedir}Part{i}/pQTL_Expsure_selection_clumping.R")
    os.chdir(Dir)
    os.system(f"sbatch sbatch_command.sh ")
    os.chdir(basedir)


##################
import os,glob
import pandas as pd

Files=glob.glob("Part*/All_significannt_CisTransexposure_BeforeQC_LDclumping.csv")
Folders=sorted([x.replace("/All_significannt_CisTransexposure_BeforeQC_LDclumping.csv","") for x in Files])
Folders=[x for x in Folders  if x!='Part1']
basedir=os.getcwd()+"/"

SubmittedFolder=['Part0', 'Part10', 'Part11', 'Part12', 'Part13', 'Part14', 'Part15', 'Part16', 'Part2', 'Part23', 'Part24', 'Part25', 'Part26', 'Part27', 'Part28', 'Part29', 'Part3', 'Part30', 'Part31', 'Part4', 'Part42', 'Part43', 'Part44', 'Part45', 'Part46', 'Part47', 'Part48', 'Part49', 'Part5', 'Part6', 'Part7', 'Part8', 'Part9']

for foder in Folders:
    Dir=f"{basedir}{foder}/"
    os.system(f"cp Decode_twosampleMR_Running.R mr_analysis_sbatch_command.sh {Dir}")
    os.chdir(Dir)
    os.system(f"sbatch mr_analysis_sbatch_command.sh ")
    os.chdir(basedir)





