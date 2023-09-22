
import pandas as pd
import numpy as np
import os,glob


rscripts=['TransExposureNoMHC_MRAnalysis_Running.R', 'TransExposure_MRAnalysis_Running.R', 'CisExposure_MRAnalysis_Running.R']

pqtltype="Biogen"

filename_dict={"daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.trios_GRCh38_UniqID.vcf.gz":"PGC3_SCZ_NoUKB",
              "ADHD2022_iPSYCH_deCODE_PGC.meta_GRCh38_UniqID.vcf.gz":"PGC_ADHD2022_iPSYCH_deCODE",
              "daner_bip_pgc3_nm_noukbiobank_GRCh38_UniqID.vcf.gz":"BIP_PGC3_noukb",
              "iPSYCH-PGC_ASD_Nov2017_GRCh38_UniqID.vcf.gz":"ASD_PGC",
              "PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf.gz":"PGC3_SCZ",
              "PGC_UKB_depression_genome-wide_GRCh38_UniqID.vcf.gz":"PGC_MDD_Depression",
              "daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01_GRCh38_UniqID.vcf.gz":"Depression_iPSYCH_2023",
              "pgcAN2.2019-07_GRCh38_UniqID.vcf.gz":"PGC_AN2"}

basedir=os.getcwd()+"/"
print(rscripts)

for rscript in rscripts:
    r_prefix=rscript.replace("_MRAnalysis_Running.R","")
    for file,name in filename_dict.items():
        file=file.replace(".gz","")
        os.system(f'''sed '/prefix="cis_exposure"/c\prefix="{pqtltype}_{name}_{r_prefix}"' {rscript} > {pqtltype}_{name}_{rscript}''')
        os.system(f'''sed -i 's/PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf/{file}/g' {pqtltype}_{name}_{rscript} ''')
        os.system(f'''sed 's/Rscript/Rscript {pqtltype}_{name}_{rscript}/' MRAnalysis_sbatch_command.sh > {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh ''')
       
        if "TransExposureNoMHCUnique" not in rscript:
            print(rscript)
            for i in range(0,50):
                Dir=f"{basedir}Part{i}/"
                os.system(f"cp {pqtltype}_{name}_{rscript} {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh {Dir}")
                os.chdir(Dir)
                os.system(f"sbatch --output={pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command-%j.out {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh ")
                os.chdir(basedir)
        os.system(f"rm {pqtltype}_{name}_{rscript} {pqtltype}_{name}_{r_prefix}_MRAnalysis_sbatch_command.sh")


