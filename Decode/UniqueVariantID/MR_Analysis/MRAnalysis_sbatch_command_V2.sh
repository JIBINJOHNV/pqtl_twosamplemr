#!/bin/bash

#SBATCH --array=0-49
#SBATCH --ntasks=4
#SBATCH --mem=260GB
#SBATCH --partition=cpu
#SBATCH --time=200:00:00

#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

source /home/jjohn1/modulefiles/anaconda3/bin/activate
module load libxml2/2.9.8-GCCcore-6.4.0
module load OpenSSL/1.1
module load R/4.1.3-foss-2021b
module load GMP/6.2.1-GCCcore-11.2.0

# Define the base directory where your folders are located
base_directory="/edgehpc/dept/human_genetics/users/jjohn1/Cis_Trans/Exposure/UKB_50KClumped/Decode/MAF_0.01/"

# Define the specific folder for this array task
folder_name="Part${SLURM_ARRAY_TASK_ID}"

# Change to the specific folder
cd "${base_directory}${folder_name}"


Rscript
