"""
Script de récupération et d'annotation taxonomique depuis la base de données NCBI Taxonomy.

Ce script lit un fichier TSV contenant une colonne de taxids, interroge la base NCBI via l'API Entrez,
et ajoute les niveaux taxonomiques (phylum, class, order, etc.) correspondants aux entrées.
Il utilise un cache local pour éviter les requêtes redondantes et permet de relancer les échecs.

Fichiers requis :
- Un fichier TSV d'entrée contenant une colonne de taxids.
- Un fichier `NCBI_api_key.yaml` avec les identifiants suivants :
  Entrez:
    email: "votre_email@domaine.com"
    api_key: "votre_clé_ncbi"

Sorties :
- Deux fichiers :
    * `*_complete_taxo.tsv` : lignes annotées avec tous les niveaux taxonomiques attendus.
    * `*_incomplete_taxo.tsv` : lignes annotées partiellement (échecs ou taxonomie incomplète).

Utilisation (en ligne de commande) :
    python get_all_taxo_from_NCBI.py --input reference_genome_summary.tsv --taxid_column taxid --batch_size 10

Auteur : PNRIA
Date : 2025-05-06
"""

import os
import pandas as pd
from joblib import Memory
from Bio import Entrez
from tqdm import tqdm
import yaml
import argparse
from wisp_light.dataset.refSeqDataset import TAXO_LEVELS

# Configuration d'Entrez
with open("NCBI_api_key.yaml", "r") as f:
    config = yaml.safe_load(f)
Entrez.email = config["Entrez"]["email"]
Entrez.api_key = config["Entrez"]["api_key"]

script_dir = os.path.dirname(os.path.abspath(__file__))  # Répertoire du script

# Dossier de cache
cache_dir = f"{script_dir}/ncbi_taxonomy_cache"
failed_taxid_file = f"{script_dir}/failed_taxids.txt"
# Initialiser Memory pour le cache
memory = Memory(cache_dir, verbose=0)

# Fonction pour récupérer les records de taxonomie
@memory.cache
def fetch_taxonomy_record_batch(taxid_list):
    """
    Récupère en batch les enregistrements taxonomiques depuis NCBI Entrez.

    Parameters
    ----------
    taxid_list : list of int
        Liste d'identifiants taxonomiques à interroger.

    Returns
    -------
    dict
        Dictionnaire associant chaque taxid à son enregistrement ou None en cas d'échec.
    """
    results = {}
    try:
        ids = ",".join(map(str, taxid_list))
        handle = Entrez.efetch(db="taxonomy", id=ids, retmode="xml")
        records = Entrez.read(handle)
        for record in records:
            tid = int(record["TaxId"])
            results[tid] = record
    except Exception as e:
        print(f"Erreur lors de la récupération du batch {taxid_list}: {e}")
        for tid in taxid_list:
            results[tid] = None
            with open(failed_taxid_file, "a") as f:
                f.write(f"{tid} taxid error (batch): {e}\n")
    return results


def retry_fail(failed_taxid_file):
    """
    Relance les requêtes pour les taxid qui ont échoué précédemment.

    Parameters
    ----------
    failed_taxid_file : str
        Chemin vers le fichier contenant les taxids échoués.

    Returns
    -------
    None
    """
    with open(failed_taxid_file) as f:
        taxid_failed = [int(line.strip()) for line in f]

    for taxid in taxid_failed:
        record = fetch_taxonomy_record_batch(taxid)
        res = 'BAD' if  not record else 'OK'
        print(f'retry for {taxid} -> ', res)


def cache_all_taxid(assembly_summary):
    """
    Met en cache les enregistrements taxonomiques de tous les taxids uniques du DataFrame.

    Parameters
    ----------
    assembly_summary : pandas.DataFrame
        DataFrame contenant une colonne 'taxid'.

    Returns
    -------
    None
    """
    all_taxid = assembly_summary['taxid']
    unique_taxid = all_taxid.unique()  # 150000 sur genome
    for taxid in tqdm(unique_taxid):
        record = fetch_taxonomy_record_batch(taxid)


def add_taxo_to_df_batch(df_in, taxid_column='taxid', batch_size=10):
    """
    Ajoute les colonnes de la lignée taxonomique à un DataFrame via des requêtes batch.

    Parameters
    ----------
    df_in : pandas.DataFrame
        DataFrame contenant une colonne de taxids.
    taxid_column : str, par défaut 'taxid'
        Nom de la colonne contenant les identifiants taxonomiques.
    batch_size : int, par défaut 10
        Nombre de taxids à traiter par batch.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame enrichi avec les colonnes de lignée taxonomique.
    failed_taxid : list of int
        Liste des taxids pour lesquels l'enregistrement a échoué.
    """

    df = df_in.copy()
    failed_taxid = []

    taxid_list = df[taxid_column].unique().tolist()
    all_records = {}

    print(f"Fetching taxonomy data in batches of {batch_size}...")

    for i in tqdm(range(0, len(taxid_list), batch_size)):
        batch = taxid_list[i:i + batch_size]
        batch_records = fetch_taxonomy_record_batch(tuple(batch))  # tuple for cache hashing
        all_records.update(batch_records)

    for index, row in tqdm(df.iterrows()):
        taxid = row[taxid_column]
        record = all_records.get(taxid)
        if record:
            lineage = record.get('LineageEx', [])
            for level in lineage[1:]:  # skip root
                rank = level['Rank']
                scientific_name = level['ScientificName']
                if rank and rank not in df.columns:
                    df[rank] = None
                df.at[index, rank] = scientific_name
        else:
            print(f"Erreur : taxid {taxid} n’a retourné aucun enregistrement.")
            failed_taxid.append(taxid)

    return df, failed_taxid

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Annoter un fichier de génomes avec la taxonomie NCBI.")
    parser.add_argument('--input', type=str, default=f'{script_dir}/reference_genome_summary.tsv',
                        help='Chemin vers le fichier TSV d’entrée contenant les taxids.')
    parser.add_argument('--taxid_column', type=str, default='taxid',
                        help='Nom de la colonne contenant les TaxIDs.')
    parser.add_argument('--batch_size', type=int, default=10,
                        help='Nombre de TaxIDs à traiter par batch.')

    args = parser.parse_args()

    df_in = pd.read_csv(args.input, sep='\t', index_col=0)
    print(df_in.head())
    df_out, failed_taxid = add_taxo_to_df_batch(df_in, taxid_column=args.taxid_column, batch_size=args.batch_size)

    print(f"Nombre total d'entrées dans df_out : {len(df_out)}")
    # Taxonomie complètes : toutes les colonnes taxonomiques présentes
    df_complete_taxo = df_out.dropna(subset=TAXO_LEVELS)

    df_incomplete_taxo = df_out[~df_out.index.isin(df_complete_taxo.index)]
    base, ext = os.path.splitext(args.input)

    print(f"Entrées avec taxonomie complète : {len(df_complete_taxo)}")
    print(f"Entrées avec taxonomie incomplète : {len(df_incomplete_taxo)}")
    df_complete_taxo.to_csv(f"{base}_complete_taxo{ext}", sep='\t', index=True)
    df_incomplete_taxo.to_csv(f"{base}_incomplete_taxo{ext}", sep='\t', index=True)



