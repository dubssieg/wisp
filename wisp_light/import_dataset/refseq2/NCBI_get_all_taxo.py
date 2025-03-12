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
        record = fetch_taxonomy_record(taxid)
        res = 'BAD' if  not record else 'OK'
        print(f'retry for {taxid} -> ', res)


def cache_all_taxid(complete_assembly_summary):
    all_taxid = complete_assembly_summary['taxid']
    unique_taxid = all_taxid.unique()  # 150000 sur genome
    for taxid in tqdm(unique_taxid):
        record = fetch_taxonomy_record(taxid)


def add_taxo_to_df(df_in, identify_name= 'taxid'):
    failed_taxid = []
    df = df_in.copy()
    for index, row in tqdm(df.iterrows()):
        taxid = row[identify_name]
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
            failed_taxid.append(taxid)
    return df, failed_taxid

if __name__ == '__main__':
    # file = f'{script_dir}/Complete_Genome_with_accession.tsv'
    # complete_assembly_summary = pd.read_csv(file, sep='\t', index_col=0)
    #
    # df, ailed_taxid  = add_taxo_to_df(complete_assembly_summary)
    # CompleteGenomeWithTaxo= pd.read_csv(f'{script_dir}/CompleteGenomeWithTaxo.csv', sep=';', index_col=0)
    df_in_file = f'{script_dir}/refseq_referent_genome_full_info.tsv'
    df_out_file = f'{script_dir}/refseq_referent_genome_with_taxo.tsv'
    df_complete_file = f'{script_dir}/complete_refseq_referent_genome_with_taxo.tsv'
    df_incomplete_file = f'{script_dir}/incomplete_refseq_referent_genome_with_taxo.tsv'

    df_in = pd.read_csv(df_in_file,  sep='\t', index_col=0)
    df_out, failed_taxid = add_taxo_to_df(df_in,identify_name='Organism Taxonomic ID' )

    df = (df_out.drop(columns=['subspecies', 'no rank', 'species group',
                               'species subgroup', 'clade', 'forma specialis', 'tribe', 'strain',
                               'suborder', 'subclass', 'subfamily', 'serotype', 'subgenus',
                               'Organism Infraspecific Names Strain',
                               'Organism Infraspecific Names Cultivar',
                               'Organism Infraspecific Names Ecotype',
                               'Organism Infraspecific Names Isolate',
                               'Organism Infraspecific Names Sex', 'Annotation BUSCO Complete ',
                               'Annotation BUSCO Single Copy ', 'Annotation BUSCO Duplicated ',
                               'Annotation BUSCO Fragmented ', 'Annotation BUSCO Missing ',
                               'Annotation BUSCO Lineage ', ]))

    df.to_csv(df_out_file, sep='\t', index=True)

    taxonomic_cols = ['phylum', 'class', 'order', 'family']
    df_complete = df.dropna(subset=taxonomic_cols)
    df_incomplete = df[~df.index.isin(df_complete.index)]
    print(len(df_complete), len(df_incomplete))
    df_complete.to_csv(df_complete_file, sep='\t', index=True)
    df_incomplete.to_csv(df_incomplete_file, sep='\t', index=True)



