#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=26
#SBATCH --time=50:00:00
#SBATCH --partition=blanca-cmbmgem
#SBATCH --output=sample-%j.out
#SBATCH --account=blanca-cmbmgem
#SBATCH --qos=blanca-cmbmgem
#SBATCH --mail-type=all
#SBATCH --mail-user=jodo9280@colorado.edu
#SBATCH --output=output.%j.out
#SBATCH --mem=15G

module purge
module load anaconda
conda activate eco_v1

python Encen_HPC.py
