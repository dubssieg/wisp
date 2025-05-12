#!/bin/bash
#SBATCH --job-name=importRefSeq	 # Nom du job
#SBATCH --chdir=/home/genouest/cnrs_umr6074/hcourtei/codes/wisp/wisp_light/build_dataset/refseq  # Répertoire de travail
#SBATCH --output=slurm/importRefSeq%j.txt     # Fichier de sortie
#SBATCH --error=slurm/importRefSeq%j.err
#SBATCH --ntasks=1               # Nombre de tâches
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00        # Durée maximale (2 jours)
#SBATCH --mem-per-cpu=20000	 # Mémoire par CPU

# sourcing conda
. /local/env/envconda.sh

# Charger l'env
. ~/envtaxo2/bin/activate

python download_refseq_from_csv.py --num_workers 5  --output_dir /projects/microtaxo/data/refseq3


