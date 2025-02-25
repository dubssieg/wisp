#!/bin/bash
#SBATCH --job-name=trainBio	 # Nom du job
#SBATCH --chdir=/home/genouest/cnrs_umr6074/hcourtei/codes/wisp/wisp_light/training  # Répertoire de travail
#SBATCH --output=trainBio.txt     # Fichier de sortie
#SBATCH --error=trainBio.err
#SBATCH --ntasks=4               # Nombre de tâches
#SBATCH --time=2-00:00:00        # Durée maximale (2 jours)
#SBATCH --mem-per-cpu=40000	 # Mémoire par CPU

# Charger l'env
. ~/envtaxo2/bin/activate

python train_val.py
