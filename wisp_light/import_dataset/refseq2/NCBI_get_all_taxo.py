import os
import pandas as pd
from joblib import Memory
from Bio import Entrez
from tqdm import tqdm


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
            f.write(f"{taxid}\n")
        return None

if __name__ == '__main__':
    print(os.getcwd())
    file = f'{script_dir}/Complete_Genome_with_accession.tsv'
    complete_assembly_summary = pd.read_csv(file, sep='\t')
    all_taxid = complete_assembly_summary['taxid']
    unique_taxid = all_taxid.unique() # 150000 sur genome
    for taxid in tqdm(unique_taxid):
        record = fetch_taxonomy_record(taxid)
