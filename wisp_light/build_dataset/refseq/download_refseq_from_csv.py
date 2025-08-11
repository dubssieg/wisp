
"""
Téléchargement et décompression de fichiers génomiques RefSeq.

Ce script lit un fichier TSV contenant des URLs FTP de génomes, télécharge
les fichiers `.fna.gz` associés, les décompresse, et met à jour un DataFrame
avec les chemins locaux et informations d'accession.

Il permet également de filtrer les génomes selon un critère donné (par exemple,
'reference genome' ou 'Complete Genome') avant le téléchargement.

Le fichier filtré est sauvegardé dans le même répertoire que le fichier CSV d'entrée,
sous le nom : <filter_by>_summary.tsv (avec les espaces remplacés par des underscores).

Utilise le multithreading pour accélérer le traitement des fichiers.

Parameters
----------
--csv_file : str
    Chemin vers le fichier TSV contenant la colonne 'ftp_path'.
    Par défaut : "assembly_summary.txt".

--output_dir : str
    Répertoire de sortie où seront placés les fichiers décompressés.
    Par défaut : "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data".

--filter_by : str
    Critère de filtrage des génomes. Doit être 'reference genome' ou 'Complete Genome'.
    Par défaut : 'reference genome'.

--num_workers : int
    Nombre de threads à utiliser pour les téléchargements/décompressions parallèles.
    Par défaut : 5.

Usage (en ligne de commande)
----------------------------
python download_refseq_from_csv.py --csv_file assembly_summary.txt \
    --output_dir /home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data \
    --filter_by reference genome \
    --num_workers 8
"""

import argparse
import os
import subprocess
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor




def update_downloaded_column(df, output_dir):
    """
    Met à jour les colonnes 'Downloaded', 'file' et 'accession' du DataFrame.

    Pour chaque ligne, vérifie si le fichier décompressé existe, extrait l'accession
    de la première ligne FASTA, et complète les colonnes.

    Parameters
    ----------
    df : pandas.DataFrame
       DataFrame contenant au moins une colonne 'ftp_path'.
    output_dir : str
       Répertoire où sont stockés les fichiers décompressés.

    Returns
    -------
    pandas.DataFrame
       Le DataFrame mis à jour avec les colonnes 'Downloaded', 'file', et 'accession'.
    """
    # Ajouter une nouvelle colonne ou mettre à jour la colonne Downloaded
    for index, row in df.iterrows():
        # Extraire l'URL FTP et créer le chemin du fichier attendu
        ftp_path = row['ftp_path']
        end_url_file = ftp_path[8:].split('/')[-1]  # On enlève 'https://' et on prend la dernière partie de l'URL
        file_path = os.path.join(output_dir, f"{end_url_file}_genomic.fna.gz")
        unzip_file_path = file_path.removesuffix(".gz")


        # Vérifier si le fichier .fna.gz existe
        if os.path.isfile(unzip_file_path):
            df.at[index, 'file'] = os.path.basename(unzip_file_path)
            df.at[index, 'Downloaded'] = True
            with open(unzip_file_path, "r", encoding='utf-8') as reader:
                # first_line = reader.readline().strip()
                first_line = reader.readline()
                accession = first_line.split('.')[0][1:]  # Extraction de l'accession
                df.at[index, 'accession']  = accession
        else:
            df.at[index, 'Downloaded'] = False
    return df


def download_and_decompress(https_path, output_dir):
    """
    Télécharge et décompresse un fichier génomique à partir d'un lien FTP.

    Utilise `wget` pour le téléchargement et `gzip` pour la décompression.
    Ignore les fichiers déjà décompressés.

    Parameters
    ----------
    https_path : str
     Lien HTTPS vers le fichier à télécharger.
    output_dir : str
     Répertoire où enregistrer les fichiers.

    Returns
    -------
    str or None
     Chemin vers le fichier décompressé, ou None si une erreur s'est produite.
    """
    ftp_path = https_path[8:]  # On enlève 'https://'
    end_url_file = ftp_path.split('/')[-1]  # Dernière partie du chemin, le nom du fichier

    file_path = os.path.join(output_dir, f"{end_url_file}_genomic.fna.gz")
    real_ftp_path = f"{ftp_path}/{end_url_file}_genomic.fna.gz"
    unzip_file_path = file_path.removesuffix(".gz")

    if os.path.exists(unzip_file_path):
        # print(f"{unzip_file_path} already exists, skipping download and decompression.")
        return unzip_file_path

    if os.path.exists(file_path):
        print(f"{file_path} seems corrupted, removing and retrying download.")
        os.remove(file_path)

    wget_command = f"wget -q -P {output_dir} {real_ftp_path}"   # Commande pour télécharger le fichier avec wget

    # Lancer wget dans un processus parallèle
    wget_process = subprocess.Popen(wget_command, shell=True)
    wget_process.wait()  # Attendre que wget ait terminé
    # Vérification si wget a rencontré une erreur
    if wget_process.returncode != 0:
        print(f"Error downloading {real_ftp_path}. Marking as failed (HTTP 404 or other issue).")
        return None  # Retourner None si échec


    decompress_command = f"gzip -f -d {file_path}"  # Commande pour décompresser le fichier téléchargé
    decompress_process = subprocess.Popen(decompress_command, shell=True)
    decompress_process.wait()  # Attendre que gzip ait terminé
    if decompress_process.returncode != 0:
        print(f"Error unzipping {real_ftp_path}. Removing corrupted file.")
        os.remove(file_path)
        return None  # Retourner None si échec

    return unzip_file_path


# Fonction pour télécharger et décompresser tous les fichiers


def download_and_decompress_all(df, output_dir, num_workers=10):
    """
    Télécharge et décompresse en parallèle les fichiers listés dans le DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame contenant les chemins FTP dans une colonne 'ftp_path'.
    output_dir : str
        Répertoire cible pour les fichiers.
    num_workers : int, optional
        Nombre de threads parallèles (default: 10).

    Returns
    -------
    list of str or None
        Liste des chemins vers les fichiers décompressés, ou None pour ceux échoués.
    """
    # Télécharger et décompresser les fichiers en parallèle
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(download_and_decompress,  row['ftp_path'], output_dir)
            for index, row in df.iterrows()]

        # Attendre que toutes les tâches soient terminées
        processed_files = []
        for future in tqdm(futures, total=len(futures), desc="Téléchargement et décompression", unit="fichier"):
            processed_files.append(future.result())

    return processed_files

def read_assembly_summary(txt_file):
    """
        Lit un fichier d'assemblage NCBI (assembly_summary.txt) et retourne un DataFrame propre.

        Ce fichier contient un en-tête commençant par le caractère `#`, que cette fonction nettoie
        pour faciliter l'accès aux colonnes.

        Parameters
        ----------
        txt_file : str or Path
            Chemin vers le fichier `assembly_summary.txt` au format tabulé (TSV).

        Returns
        -------
        pandas.DataFrame
            Un DataFrame contenant les métadonnées d'assemblage, avec le nom de la première
            colonne (`#assembly_accession`) nettoyé (le `#` est retiré).
        """
    completed_df = pd.read_csv(txt_file, sep="\t", header=1)
    completed_df.columns = [col.lstrip('#') for col in completed_df.columns]
    print(f"Lecture en Dataframe de assembly_summary.txt. Contient {len(completed_df)} génomes ")

    return completed_df
def filter_assembly_summary_df(df, filter_by):
    if filter_by == "reference genome":
        by_column = 'refseq_category'
    elif filter_by == "Complete Genome":
        by_column = 'assembly_level'
    else:
        raise NotImplemented (f"Filtrage par {filter_by} non implémenté")
    filter_df = df[df[by_column]==filter_by].reset_index(drop=True)
    print(f"Filtrage par '{by_column}' = '{filter_by}' : {len(filter_df)} génomes retenus.")

    return filter_df

if __name__ == "__main__":
    #
    parser = argparse.ArgumentParser(description="Téléchargement et décompression de fichiers génomiques.")
    parser.add_argument("--output_dir", type=str,
                        help="Répertoire de sortie pour les fichiers téléchargés et décompressés",
                        default='/home/hcourtei/Projects/MicroTaxo/codes/data/refseq_data')
    parser.add_argument("--csv_file", type=str, help="Fichier CSV contenant les chemins FTP",
                        default="assembly_summary.txt")
    parser.add_argument("--filter_by", type=str, help="filtrage des genomes: 'reference genome' ou 'Complete Genome'",
                        default="reference genome")
    parser.add_argument("--num_workers", type=int, default=5, help="Nombre de workers pour le téléchargement et la décompression (par défaut 5)")

    args = parser.parse_args()
    # args.output_dir = '/wisp/wisp_light/import_dataset/refseq/out_refseq'
    # completed_csv_out = "assembly_summary_Complete_Genome.csv"
    completed_df = read_assembly_summary(args.csv_file)
    filter_df = filter_assembly_summary_df(completed_df, filter_by=args.filter_by)
    filter_name = args.filter_by.replace(" ", "_")
    filtered_csv_path = os.path.join( os.path.dirname(args.csv_file), f"{filter_name}_summary.tsv")
    filter_df.to_csv(filtered_csv_path, sep="\t", index=False)
    print(f"Fichier filtré sauvegardé : {filtered_csv_path}")
    # # Créer le répertoire de sortie si nécessaire
    os.makedirs(args.output_dir, exist_ok=True)
    #
    # # Lancer le téléchargement et la décompression
    # processed_files = download_and_decompress_all(completed_df, args.output_dir, num_workers=args.num_workers)
    # # #
    # print(f"{len(processed_files)} fichiers traités.")



''' Batch_request

import requests
import yaml
from Bio import SeqIO, Entrez
from io import StringIO


# batch = accessions[i:i + batch_size]
# with Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text") as taxo_handle:
#     # records = SeqIO.read(taxo_handle, 'genbank')
#     records = SeqIO.parse(taxo_handle, 'genbank')
#     for idx , record in enumerate(records):
#         accession = record.id.split('.')[0]
#         taxonomy = record.annotations.get('taxonomy', [])
#         organism = record.annotations.get('organism', "Unknown Organism")

'''