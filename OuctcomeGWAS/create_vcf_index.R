#.libPaths("/home/jjohn1/modulefiles/R4.3.0-foss-2021b")
#library(argparse)

.libPaths("/edgehpc/dept/human_genetics/jrobins_legacy/r_libs_4.1.0")
#.libPaths("/home/jjohn1/modulefiles/R4.3.0-foss-2021b")

suppressWarnings(suppressPackageStartupMessages({
    library(glue)
    library(gwasvcf)
    library(VariantAnnotation)
    library(dplyr)
    library(magrittr)
    library(VariantAnnotation)
    library(TwoSampleMR)
    library(MRInstruments)
    library(gwasglue)
    library(data.table)
    library(ieugwasr)
    library(plyr)
    library(mrpipeline)
    library(data.table)
}))


outcome_gwassumstat <- c("ADHD2022_iPSYCH_deCODE_PGC.meta_GRCh38_UniqID.vcf.gz","iPSYCH-PGC_ASD_Nov2017_GRCh38_UniqID.vcf.gz",
                         "PGC_UKB_depression_genome-wide_GRCh38_UniqID.vcf.gz","daner_bip_pgc3_nm_noukbiobank_GRCh38_UniqID.vcf.gz",
                         "PGC3_SCZ_wave3.european.autosome.public.v3_GRCh38_UniqID.vcf.gz",
                        "daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01_GRCh38_UniqID.vcf.gz",
                        "pgcAN2.2019-07.vcf._GRCh38_UniqID.vcf.gz")

path <- "/mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_Processed_Outcome_GWAS_UniqIds/"

for (outcomevcf in outcome_gwassumstat) {
      full_path_vcf <- glue(path, outcomevcf)
      outcomevcf_index <- glue("{sub('.gz$', '', full_path_vcf)}.index")
      create_rsidx_index_from_vcf(full_path_vcf, outcomevcf_index)
      }

#srun -n 12 --mem=30G --time=1:01:00  --pty bash
#module load Python/3.8.8-GCCcore-11.2.0
#module load libxml2/2.9.8-GCCcore-6.4.0
#module load OpenSSL/1.1
#source /edgehpc/dept/human_genetics/users/jjohn1/Discovery_Genewise_Merged/gwas2vcf/env/bin/activate
#export PATH="/edgehpc/dept/human_genetics/users/jjohn1/Discovery_Genewise_Merged/gwas2vcf/env/bin:$PATH"
