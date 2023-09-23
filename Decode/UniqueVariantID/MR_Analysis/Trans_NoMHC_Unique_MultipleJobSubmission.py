import pandas as pd
import numpy as np
import os,glob

##For this clum results from the 50 folder neeeds to merge and create Biogen_All_significannt_trans_exposure_AfterQC_LDclumping_NoMHC_Unique.csv
pqtltype="Decode"

rscripts=glob.glob("TransExposureNoMHCUnique_MRAnalysis_Running.R")
basedir=os.getcwd()+"/"
print(rscripts)

file=f"{basedir}Decode_CombinedResultsfromAllBatches/{pqtltype}_All_significannt_trans_exposure_AfterQC_LDclumping_NoMHC_Unique.csv"
df=pd.read_csv(file)
df=df.sort_values(by=["exposure","ID"]).reset_index(drop="index")
exposure_list=list(df["exposure"].unique())

filename_dict={"daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.trios_GRCh38_UniqID.vcf.gz":"PGC3_SCZ_NoUKB",
              "ADHD2022_iPSYCH_deCODE_PGC.meta_GRCh38_UniqID.vcf.gz":"PGC_ADHD2022_iPSYCH_deCODE",
              "daner_bip_pgc3_nm_noukbiobank_GRCh38_UniqID.vcf.gz":"BIP_PGC3_noukb",
              "iPSYCH-PGC_ASD_Nov2017_GRCh38_UniqID.vcf.gz":"ASD_PGC",
              "PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf.gz":"PGC3_SCZ",
              "daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01_GRCh38_UniqID.vcf.gz":"Depression_iPSYCH_2023",
              "pgcAN2.2019-07_GRCh38_UniqID.vcf.gz":"PGC_AN2"}

result_df=pd.DataFrame()

for rscript in rscripts:
    r_prefix=rscript.replace("_MRAnalysis_Running.R","")
    
    for file,name in filename_dict.items():
        file=file.replace(".gz","")
        os.system(f'''sed '/prefix="cis_exposure"/c\prefix="{pqtltype}_{name}_{r_prefix}"' {rscript} > {pqtltype}_{name}_{rscript}''')
        os.system(f'''sed -i 's/PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf/{file}/g' {pqtltype}_{name}_{rscript} ''')
        os.system(f'''sed 's/Rscript/Rscript {pqtltype}_{name}_{rscript}/' MRAnalysis_sbatch_command.sh > {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh ''')

        n=0
        for i in range(0,50):
            temp_df=df[df["exposure"].isin(exposure_list[n:n+41])]
            result_df=pd.concat([temp_df,result_df])
            n=n+41
            
            Dir=f"{basedir}Part{i}/"
            temp_df.to_csv(f"{Dir}All_significannt_trans_exposure_AfterQC_LDclumping_NoMHC_Unique.csv")
            os.system(f"cp {pqtltype}_{name}_{rscript} {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh {Dir}")
            os.chdir(Dir)
            os.system(f"sbatch --output={pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command-%j.out {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh ")
            os.chdir(basedir)
        os.system(f"rm {pqtltype}_{name}_{rscript} {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh")

