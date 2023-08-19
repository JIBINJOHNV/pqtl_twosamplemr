#!/bin/bash

#SBATCH --ntasks=4
#SBATCH --time=00:01:00
#SBATCH --mem=120GB
#SBATCH --partition=cpu
#SBATCH --time=100:00:00


source /home/jjohn1/modulefiles/anaconda3/bin/activate
module load SAMtools/1.15-GCC-11.2.0
module load PLINK/1.9b_6.21-x86_64
module load R/4.1.3-foss-2021b





Rscript Decode_pQTL_Expsure_selection_clumping_UniqVariantID.R
