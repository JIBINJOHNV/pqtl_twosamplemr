import pandas as pd
import numpy as np
import os


# Load your original DataFrame
df = pd.read_csv("/edgehpc/dept/human_genetics/users/jjohn1/Cis_Trans/Exposure/UKB_50KClumped/Biogen/PGC3_SCZ/Unique_variant_Id/All_Biogen_Exposure_TSS_chr_index_3k.csv")


# Split the DataFrame into almost equal-sized parts
split_dfs = np.array_split(df, 50)

basedir=os.getcwd()+"/"

# Save each part to separate files
for i, split_df in enumerate(split_dfs):
    os.system(f"mkdir -p Part{i}")
    Dir=f"{basedir}Part{i}/"
    split_df.to_csv(f"{Dir}All_Biogen_Exposure_TSS_chr_index_3k_Part{i}.csv", index=False)
    os.system(f"cp Biogen_sbatch_command_clumping.sh Biogen_pQTL_Expsure_selection_clumping.R {Dir}")
    os.system(f'''sed -e "s|All_Biogen_Exposure_TSS_chr_index_3k_Part0.csv|{Dir}All_Biogen_Exposure_TSS_chr_index_3k_Part{i}.csv|" Biogen_pQTL_Expsure_selection_clumping.R > {basedir}Part{i}/Biogen_pQTL_Expsure_selection_clumping.R''')
    os.chdir(Dir)
    os.system(f"sbatch Biogen_sbatch_command_clumping.sh ")
    os.chdir(basedir)
