import os
import pandas as pd
from joblib import Memory
from Bio import Entrez
from tqdm import tqdm
import json


# Configuration d'Entrez
Entrez.email = "hermann.courteille@inria.fr"
Entrez.api_key = "b55513ab1634ec527ccf1ec084f3b1c78108"
script_dir = os.path.dirname(os.path.abspath(__file__))  # Répertoire du script

# Dossier de cache
cache_dir = f"{script_dir}/ncbi_taxonomy_cache"
failed_taxid_file = f"{script_dir}/failed_taxids.txt"
# Initialiser Memory pour le cache
memory = Memory(cache_dir, verbose=0)

# Fonction pour récupérer les records de taxonomie
@memory.cache
def fetch_taxonomy_record(taxid):
    try:
        handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
        records = Entrez.read(handle)
        # print(f"TaxID {taxid} ajouté au cache")
        return records
    except Exception as e:
        print(f"Erreur lors de la récupération du TaxID {taxid}: {e}")
        with open(failed_taxid_file, "a") as f:
            f.write(f"{taxid} taxid error : {e}\n")
        return None

def retry_fail(failed_taxid_file):
    with open(failed_taxid_file) as f:
        taxid_failed = [int(line.strip()) for line in f]

    for taxid in taxid_failed:
        print('retry for : ', taxid)
        record = fetch_taxonomy_record(taxid)

def cache_all_taxid(complete_assembly_summary):
    all_taxid = complete_assembly_summary['taxid']
    unique_taxid = all_taxid.unique()  # 150000 sur genome
    for taxid in tqdm(unique_taxid):
        record = fetch_taxonomy_record(taxid)


def add_taxo_to_df(complete_assembly_summary, sep=';'):
    failed_taxid = []
    df = complete_assembly_summary.copy()
    for index, row in tqdm(df.iterrows()):
        taxid = row['taxid']
        record = fetch_taxonomy_record(taxid)
        if record:
            lineage = record[0]['LineageEx']
            for level in lineage[1:]:
                rank = level['Rank']
                scientific_name = level['ScientificName']
                if rank and rank not in df.columns:
                    df[rank] = None

                df.at[index, rank] = scientific_name
        else:
            print(f"error taxid {taxid} empty record {record}")
            failed_taxid.append((taxid))
    df.to_csv(f'{script_dir}/CompleteGenomeWithTaxo.csv',sep=sep, index=True)
    return df, failed_taxid

if __name__ == '__main__':
    file = f'{script_dir}/Complete_Genome_with_accession.tsv'
    complete_assembly_summary = pd.read_csv(file, sep='\t', index_col=0)

    # df, ailed_taxid  = add_taxo_to_df(complete_assembly_summary)
    CompleteGenomeWithTaxo= pd.read_csv(f'{script_dir}/CompleteGenomeWithTaxo.csv', sep=';', index_col=0)

    columns_to_drop = [
        'suborder', 'subclass', 'clade', 'species', 'no rank', 'species group',
        'species subgroup', 'subfamily', 'subgenus', 'subspecies', 'serogroup',
        'tribe', 'biotype', 'forma specialis', 'serotype', 'strain', 'pathogroup', 'ftp_path', 'species_taxid', 'file'
    ]
    missing_values = CompleteGenomeWithTaxo.isna().sum()
    print("missing_values", missing_values)
    df = CompleteGenomeWithTaxo.drop(columns=columns_to_drop)

    df_cleaned = df.dropna(subset=['kingdom', 'phylum', 'class', 'order'])
    df_no_family = df_cleaned[df_cleaned['family'].isna()]
    print(df_no_family.to_markdown())
    df_with_family = df_cleaned[df_cleaned['family'].notna()]
    print(f" no family len {len(df_no_family)} with {len(df_with_family)}")

