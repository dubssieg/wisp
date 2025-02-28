#!/bin/bash
#SBATCH --job-name=trainBio	 # Nom du job
#SBATCH --chdir=/home/genouest/cnrs_umr6074/hcourtei/codes/wisp/wisp_light/training  # Répertoire de travail
#SBATCH --output=slurm/trainBio_%j.txt     # Fichier de sortie
#SBATCH --error=slurm/trainBio_%j.err
#SBATCH --ntasks=1               # Nombre de tâches
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00        # Durée maximale (2 jours)
#SBATCH --mem-per-cpu=20000	 # Mémoire par CPU

##SBATCH --mem=80G
##SBATCH --mail-user=siegfried.dubois@irisa.fr
##SBATCH --mail-type=END,FAIL

# sourcing conda
. /local/env/envconda.sh

# variables
#WD=$(pwd)
#ENV=$WD"/wisp-env"
# Charger l'env
. ~/envtaxo2/bin/activate

python train_val.py

echo "Le numéro du job SLURM est : $SLURM_JOB_ID"

scontrol show job $SLURM_JOB_ID > slurm/job_config_${SLURM_JOB_ID}.txt


