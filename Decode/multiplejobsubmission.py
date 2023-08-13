
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


#################### Trans Exposure MR NNoMHC analysis multiple job
import pandas as pd
import numpy as np
import os
basedir=os.getcwd()+"/"
# Save each part to separate files
for i in range(0,50):
    Dir=f"{basedir}Part{i}/"
    os.system(f"cp Decode_TransExposure_NoMHC_twosampleMR_Running.R TransExposure_NoMHC_mr_analysis_sbatch_command.sh {Dir}")
    os.chdir(Dir)
    os.system(f"sbatch TransExposure_NoMHC_mr_analysis_sbatch_command.sh ")
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



################## Meerge results of cis_exposure_ MR analysis results
import os,glob
import pandas as pd

unique_files=glob.glob("Part0/cis_exposure_*")
unique_files=[x.split("/")[1] for x in unique_files ]

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





