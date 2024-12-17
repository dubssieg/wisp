#!/bin/bash

# Vérifier que group_id est fourni comme argument
if [ -z "$1" ]; then
    echo "Usage: bash submit_get_taxo_by_batch.sh <group_id>"
    exit 1
fi

# Vérifier si group_id est un entier
group_id=$1
if ! [[ "$group_id" =~ ^[0-9]+$ ]]; then
    echo "Error: group_id must be an integer."
    exit 1
fi

# Soumettre le job avec SLURM
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=GetTaxo       # Nom du job
#SBATCH --chdir=/home/genouest/cnrs_umr6074/hcourtei/codes/wisp/scripts  # Répertoire de travail
#SBATCH --output=GetTaxo.txt     # Fichier de sortie
#SBATCH --ntasks=2               # Nombre de tâches
#SBATCH --time=2-00:00:00        # Durée maximale (2 jours)
#SBATCH --mem-per-cpu=2000       # Mémoire par CPU

# Charger Python si nécessaire
. ~/envtaxo2/bin/activate

# Exécuter le script Python avec le paramètre group_id
python get_taxo_by_batch.py $group_id
EOT
