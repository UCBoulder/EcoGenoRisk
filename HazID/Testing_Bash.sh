#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=100:00:00
#SBATCH --partition=blanca-cmbmgem
#SBATCH --output=sample-%j.out
#SBATCH --account=blanca-cmbmgem
#SBATCH --qos=blanca-cmbmgem
#SBATCH --mail-type=all
#SBATCH --mail-user=jodo9280@colorado.edu
#SBATCH --output=output.%j.out

module purge
module load anaconda
conda activate eco_v1

python _1_Library_Calling_Script.py
