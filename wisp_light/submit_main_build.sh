#!/bin/bash

# Soumettre le job avec SLURM
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=Build_refseq	 # Nom du job
#SBATCH --chdir=/home/genouest/cnrs_umr6074/hcourtei/codes/wisp/workspace  # Répertoire de travail
#SBATCH --output=Build_refseq.txt     # Fichier de sortie
#SBATCH --error=Build_refseq.err
#SBATCH --ntasks=2               # Nombre de tâches
#SBATCH --time=2-00:00:00        # Durée maximale (2 jours)
#SBATCH --mem-per-cpu=40000	 # Mémoire par CPU

# Charger Python si nécessaire
. ~/envtaxo2/bin/activate

# Exécuter le script Python avec le paramètre group_id
python main.py build refseq /groups/microtaxo/data/refseq_with_taxo/
EOT
