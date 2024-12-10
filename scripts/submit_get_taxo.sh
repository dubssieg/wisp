#!/bin/bash
ID_GROUP=$1
#SBATCH --job-name=GetTaxo_"$ID_GROUP"      # Nom du job
#SBATCH --chdir=/home/genouest/cnrs_umr6074/hcourtei/codes/wisp/scripts     # Répertoire de travail
#SBATCH --output=GetTaxo_"$ID_GROUP".txt # Fichier de sortie
#SBATCH --ntasks=2                   # Nombre de tâches
#SBATCH --time=2-00:00:00            # Durée maximale (2 jours)
#SBATCH --mem-per-cpu=2000          # Mémoire par CPU



# Vérification que le paramètre est un entier
if ! [[ "$ID_GROUP" =~ ^[0-9]+$ ]]; then
  echo "Erreur : le paramètre fourni n'est pas un entier valide."
  exit 1
fi

# Lancer le script avec le paramètre entier comme argument
bash my_script.sh "$ID_GROUP"

