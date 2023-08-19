#!/bin/bash

#SBATCH --ntasks=4
#SBATCH --mem=260GB
#SBATCH --partition=cpu
#SBATCH --time=200:00:00


source /home/jjohn1/modulefiles/anaconda3/bin/activate
module load libxml2/2.9.8-GCCcore-6.4.0
module load OpenSSL/1.1
module load R/4.1.3-foss-2021b
module load GMP/6.2.1-GCCcore-11.2.0



Rscript 
