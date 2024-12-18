import time

from natsort import natsorted
import logging
from tqdm import tqdm
from argparse import ArgumentParser
import os
import pandas as pd
import matplotlib.pyplot as plt




def scan_for_taxo(datadir):
    # Initialiser une liste pour stocker les informations extraites
    data = []
    subgroups = natsorted(os.listdir(datadir))
    print("GROUP in this data dir" , subgroups)
    # Parcourir tous les sous-répertoires et fichiers
    for root, dirs, files in tqdm(os.walk(datadir)):
        for file_name in files:
            if file_name.endswith(".fna") and not file_name.startswith("GCF_"):
                # Vérifier si le fichier suit le format taxonomique attendu
                parts = file_name.split("_")
                if len(parts) >= 6:  # Il doit y avoir au moins 7 éléments dans le nom
                    regne = parts[0]
                    phylum = parts[1]
                    ordre = parts[2]
                    famille = parts[3]
                    genre = parts[4]
                    espece_with_extension = parts[5]

                    # Identifier l'ID (s'il existe)
                    try:
                        id_fichier = int(parts[6].replace(".fna", "")) if len(parts) > 6 else 0
                    except ValueError:
                        id_fichier = 0

                    # Supprimer l'extension `.fna` pour l'espèce si ID est 0
                    if id_fichier == 0:
                        espece = espece_with_extension.replace(".fna", "")
                    else:
                        espece = espece_with_extension
                        # Ajouter les informations dans la liste
                    group_name = os.path.basename(root)  # Nom du répertoire contenant le fichier
                    data.append([group_name, regne, phylum, ordre, famille, genre, espece, id_fichier])


    # Créer un DataFrame pandas à partir des données
    columns = ["Groupe", "Règne", "Phylum", "Ordre", "Famille", "Genre", "Espèce", "ID"]
    print("Extraction terminée pour les colonnes ", ' '.join(columns))

    raw_df = pd.DataFrame(data, columns=columns)
    return raw_df


def compute_counts(raw_df):
    # Étape 1 : Compter les occurrences totales par combinaison
    counts = (
        raw_df.groupby(["Règne", "Phylum", "Ordre", "Famille", "Genre", "Espèce", "Groupe"])
        .size()
        .reset_index(name="Effectif")  # Ajout de la colonne Effectif
    )
    # Étape 2 : Additionner les effectifs sur tous les groupes
    total_counts_over_group = (
        counts.groupby(["Règne", "Phylum", "Ordre", "Famille", "Genre", "Espèce"])["Effectif"]
        .sum()
        .reset_index()
    )
    family_counts = (
        counts.groupby(["Règne", "Phylum", "Ordre", "Famille"])["Effectif"]
        .sum()
        .reset_index()
    )
    return counts, total_counts_over_group, family_counts

def plot_family_counts(family_counts):
    plt.figure(figsize=(10, 6))
    plt.bar(family_counts["Famille"], family_counts["Effectif"], color="skyblue")
    plt.xticks(rotation=90)  # Faire tourner les labels de l'axe des x pour plus de lisibilité
    plt.xlabel('Famille')
    plt.ylabel('Effectif total')
    plt.title('Effectif total par Famille')
    plt.tight_layout()  # Pour éviter que le texte ne soit coupé
    plt.show()


if __name__ == '__main__':

    datadir = "/groups/microtaxo/data/refseq_unzip_with_taxo" # "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_unzip_with_taxo"  #    # #
    stat_dir = os.path.join(os.path.dirname(datadir), 'stat_'+os.path.basename(datadir))
    os.makedirs(stat_dir, exist_ok=True)
    raw_df = scan_for_taxo(datadir)

    counts, total_counts_over_group, family_counts = compute_counts(raw_df)
    raw_df.to_csv(os.path.join(stat_dir,"raw_taxonomy_data.csv"), index=False, sep= ";")
    total_counts_over_group.to_csv(os.path.join(stat_dir,"total_counts_over_group.csv"), index=False, sep=";")
    family_counts.to_csv(os.path.join(stat_dir, "family_counts.csv"), index=False, sep=";")
    print(f"All stat save in {stat_dir}")

    plot_family_counts(family_counts)