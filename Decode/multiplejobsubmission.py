
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

