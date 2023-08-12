
###Select Exposure and clumping
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


#################### Cis Exposure MR analysis multiple job

import pandas as pd
import numpy as np
import os

basedir=os.getcwd()+"/"
# Save each part to separate files
for i in range(0,50):
    Dir=f"{basedir}Part{i}/"
    os.system(f"cp mr_analysis_sbatch_command.sh  Decode_CisExposure_twosampleMR_Running.R {Dir}")
    os.chdir(Dir)
    os.system(f"sbatch mr_analysis_sbatch_command.sh ")
    os.chdir(basedir)

#################### Trans Exposure MR analysis multiple job

import pandas as pd
import numpy as np
import os
basedir=os.getcwd()+"/"
# Save each part to separate files
for i in range(0,50):
    Dir=f"{basedir}Part{i}/"
    os.system(f"cp Decode_TransExposure_twosampleMR_Running.R TransExposure_mr_analysis_sbatch_command.sh {Dir}")
    os.chdir(Dir)
    os.system(f"sbatch TransExposure_mr_analysis_sbatch_command.sh ")
    os.chdir(basedir)





#################### Trans NoMHC Unique Exposure MR analysis multiple job
import pandas as pd
import numpy as np
import os

basedir=os.getcwd()+"/"

# Save each part to separate files
for i in range(0,50):
    Dir=f"{basedir}Part{i}/"
    os.system(f"cp Decode_TransExposure_NoMHC_Unique_twosampleMR_Running.R TransExposure_NoMHC_Unique_mr_analysis_sbatch_command.sh {Dir}")
    os.chdir(Dir)
    df=pd.read_csv("All_significannt_trans_exposure_AfterQC_LDclumping_MHCRemoval.csv")
    count_df=df.groupby("ID")["seqnames"].count().reset_index()
    count_uniq=count_df[count_df["seqnames"]==1]
    count_uniq=count_uniq.rename(columns={'seqnames':"Count"})
    df=pd.merge(df,count_uniq,on="ID").drop("Count",axis=1)
    df.to_csv("All_significannt_trans_exposure_AfterQC_LDclumping_NoMHC_Unique.csv")
    os.system(f"sbatch TransExposure_NoMHC_Unique_mr_analysis_sbatch_command.sh ")
    os.chdir(basedir)





################## Meerge results from above job
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





