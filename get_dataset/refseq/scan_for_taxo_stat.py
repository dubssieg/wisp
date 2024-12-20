import os
import argparse
import pandas as pd
from natsort import natsorted
from tqdm import tqdm


def scan_for_taxo(datadir):
    # Initialiser une liste pour stocker les informations extraites
    data = []
    subgroups = natsorted(os.listdir(datadir))
    print(f"GROUP in this data dir {subgroups}" )
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
    # from 'root', 'domain', 'phylum', 'group', 'order', 'family'
    print(f"Extraction terminée pour les colonnes {' '.join(columns)}", )

    raw_df = pd.DataFrame(data, columns=columns)
    return raw_df



if __name__ == '__main__':

    # parser = argparse.ArgumentParser(description="Taxonomy Data Processor")
    # parser.add_argument("--datadir", type=str, help="Path to the data directory",
    #                     default="/groups/microtaxo/data/refseq_with_taxo")
    # args = parser.parse_args()

    datadir = "/groups/microtaxo/data/refseq_with_taxo" # "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_unzip_with_taxo"  #    # #
    raw_df = scan_for_taxo(datadir)
    output_stat_file = os.path.join(os.path.dirname(datadir),"raw_taxonomy_refseq_data.csv")
    raw_df.to_csv(output_stat_file, index=False, sep= ";")
    print(f"All stat save in {output_stat_file}")