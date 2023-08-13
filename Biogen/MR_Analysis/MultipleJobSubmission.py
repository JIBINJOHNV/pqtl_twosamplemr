#################### Cis Exposure MR analysis multiple job
import pandas as pd
import numpy as np
import os

basedir=os.getcwd()+"/"
# Save each part to separate files
for i in range(0,50):
    Dir=f"{basedir}Part{i}/"
    os.system(f"cp CisExposure_mr_analysis_sbatch_command.sh  Biogen_CisExposure_twosampleMR_Running.R {Dir}")
    os.chdir(Dir)
    os.system(f"sbatch CisExposure_mr_analysis_sbatch_command.sh ")
    os.chdir(basedir)

#################### Trans Exposure MR analysis multiple job
import pandas as pd
import numpy as np
import os
basedir=os.getcwd()+"/"
# Save each part to separate files
for i in range(0,50):
    Dir=f"{basedir}Part{i}/"
    os.system(f"cp Biogen_TransExposure_twosampleMR_Running.R TransExposure_mr_analysis_sbatch_command.sh {Dir}")
    os.chdir(Dir)
    os.system(f"sbatch TransExposure_mr_analysis_sbatch_command.sh ")
    os.chdir(basedir)

#################### Trans Exposure MR NoMHC analysis multiple job
import pandas as pd
import numpy as np
import os
basedir=os.getcwd()+"/"
# Save each part to separate files
for i in range(0,50):
    Dir=f"{basedir}Part{i}/"
    os.system(f"cp Biogen_TransExposure_NoMHC_twosampleMR_Running.R TransExposure_NoMHC_mr_analysis_sbatch_command.sh {Dir}")
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
    os.system(f"cp TransExposure_NoMHC_Unique_mr_analysis_sbatch_command.sh Biogen_TransExposure_NoMHC_Unique_twosampleMR_Running.R {Dir}")
    os.chdir(Dir)
    df=pd.read_csv("All_significannt_trans_exposure_AfterQC_LDclumping_MHCRemoval.csv")
    count_df=df.groupby("ID")["seqnames"].count().reset_index()
    count_uniq=count_df[count_df["seqnames"]==1]
    count_uniq=count_uniq.rename(columns={'seqnames':"Count"})
    df=pd.merge(df,count_uniq,on="ID").drop("Count",axis=1)
    df.to_csv("All_significannt_trans_exposure_AfterQC_LDclumping_NoMHC_Unique.csv")
    os.system(f"sbatch TransExposure_NoMHC_Unique_mr_analysis_sbatch_command.sh ")
    os.chdir(basedir)
