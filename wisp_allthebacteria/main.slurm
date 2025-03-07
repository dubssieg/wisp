#!/bin/bash
#SBATCH --job-name=clx_atb
#SBATCH --chdir=/home/genouest/cnrs_umr6074/cleroux2/microtaxo/wisp/wisp_allthebacteria
#SBATCH --output=slurm/atb_%j.txt     # Fichier de sortie
#SBATCH --error=slurm/atb_%j.err
#SBATCH --ntasks=1               # Nombre de tâches
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00        # Durée maximale (2 jours)
#SBATCH --mem-per-cpu=20000	 # Mémoire par CPU
##SBATCH --mem=40G
##SBATCH --mail-user=cyrille.leroux@irisa.fr
##SBATCH --mail-type=END,FAIL


source ~/.bashrc
conda activate wisp

python main.py --create-db
