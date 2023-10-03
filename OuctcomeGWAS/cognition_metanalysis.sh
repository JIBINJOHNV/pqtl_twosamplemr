/mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/Cognition/generic-metal/METAL/build/bin/metal

#THIS SCRIPT EXECUTES AN ANALYSIS OF EIGHT STUDIES
#THE RESULTS FOR EACH STUDY ARE STORED IN FILES Inputfile1.txt THROUGH Inputfile8.txt

#LOAD THE FIRST EIGHT INPUT FILES

# UNCOMMENT THE NEXT LINE TO ENABLE GenomicControl CORRECTION
# GENOMICCONTROL ON

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER ID
ALLELE REF ALT
EFFECT ES
PVALUE LP 
WEIGHT SS
SEPARATOR  TAB
AVERAGEFREQ ON
MINMAXFREQ ON
PROCESS /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/Cognition/Cognition_UKB/regenie_ukb_step2_linear_model2_fiall_baseline.regenie.summary_GRCh38_UniqID_ForMeta.tsv

# === THE SECOND INPUT FILE HAS THE SAME FORMAT AND CAN BE PROCESSED IMMEDIATELY ===
PROCESS /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/Cognition/Cognition_Savage_iq/Cognition_Savage_iq_metaanalysis_noUKB_GRCh38_UniqID_ForMeta.tsv 


OUTFILE Cognition_Meta_GWAS_without_UKBPP .tbl
ANALYZE HETEROGENEITY

QUIT
