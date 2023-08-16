
import pandas as pd
import numpy as np
import os,glob


rscripts=glob.glob("*_MRAnalysis_Running.R")
pqtltype="Biogen"

filename_dict={"ADHD2022_iPSYCH_deCODE_PGC.meta_GRCh38_UniqID.vcf.gz":"PGC_ADHD2022_iPSYCH_deCODE",
              "daner_bip_pgc3_nm_noukbiobank_GRCh38_UniqID.vcf.gz":"BIP_PGC3_noukb",
              "iPSYCH-PGC_ASD_Nov2017_GRCh38_UniqID.vcf.gz":"ASD_PGC",
              "PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf.gz":"PGC3_SCZ",
              "PGC_UKB_depression_genome-wide_GRCh38_UniqID.vcf.gz":"PGC_MDD_Depression"}basedir=os.getcwd()+"/"

for rscript in rscripts:
    r_prefix=rscript.replace("_MRAnalysis_Running.R","")
    for file,name in filename_dict.items():
        file=file.replace(".gz","")
        os.system(f'''sed '/prefix="cis_exposure"/c\prefix="{pqtltype}_{name}_{r_prefix}"' {rscript} > {pqtltype}_{name}_{rscript}''')
        os.system(f'''sed -i 's/PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf/{file}/g' {pqtltype}_{name}_{rscript} ''')
        os.system(f'''sed 's/Rscript/Rscript {pqtltype}_{name}_{rscript}/' MRAnalysis_sbatch_command.sh > {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh ''')
        
        if "TransExposureNoMHCUnique" in rscript:
            for i in range(0,50):
                Dir=f"{basedir}Part{i}/"
                os.system(f"cp {pqtltype}_{name}_{rscript} {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh {Dir}")
                os.chdir(Dir)
                df=pd.read_csv("All_significannt_trans_exposure_AfterQC_LDclumping_MHCRemoval.csv")
                count_df=df.groupby("ID")["seqnames"].count().reset_index()
                count_uniq=count_df[count_df["seqnames"]==1]
                count_uniq=count_uniq.rename(columns={'seqnames':"Count"})
                df=pd.merge(df,count_uniq,on="ID").drop("Count",axis=1)
                df.to_csv("All_significannt_trans_exposure_AfterQC_LDclumping_NoMHC_Unique.csv")
                os.system(f"sbatch --output={pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command-%j.out {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh  ")
                os.chdir(basedir)
        else:
            for i in range(0,50):
                Dir=f"{basedir}Part{i}/"
                os.system(f"cp {pqtltype}_{name}_{rscript} {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh {Dir}")
                os.chdir(Dir)
                os.system(f"sbatch --output={pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command-%j.out {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh ")
                os.chdir(basedir)

