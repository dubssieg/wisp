import os
import gzip
import shutil
import time
from natsort import natsorted
from Bio import Entrez, SeqIO
from tqdm import tqdm
import logging

# Configuration du logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from Bio import Entrez, SeqIO
import os
import gzip
import shutil
from tqdm import tqdm


def get_taxo_by_batch(input_dir: str, subgroup: str, batch_size: int = 100) -> None:
    """
    Get taxonomy name for genomes from NCBI using batch requests.

    Args:
        input_dir (str): Path to folder with raw genome files (.fna.gz).
        subgroup (str): Subgroup name for processing specific folder.
        batch_size (int): Number of accessions to fetch in a single request.
    """
    # Configurer les paramètres Bio.Entrez
    Entrez.email = "hermann.courteille@inria.fr"  # Remplacez par votre email
    Entrez.api_key = "b55513ab1634ec527ccf1ec084f3b1c78108"
    Entrez.max_tries = 5
    Entrez.sleep_between_tries = 15
    nb_convert = 0

    output_dir = f"{input_dir}_unzip_with_taxo"
    input_dir_sub = os.path.join(input_dir, subgroup)
    genome_file_list = [f for f in os.listdir(input_dir_sub) if f.endswith('.gz')]

    accession_map = {}
    decompressed_files = []
    # Étape 1 : Décompression des fichiers et collecte des accessions
    for genome_file in tqdm(genome_file_list, desc="Decompressing files"):
        input_filepath = os.path.join(input_dir_sub, genome_file)
        decompressed_path = os.path.join(output_dir, subgroup, genome_file[:-3])  # Retirer .gz

        try:
            os.makedirs(os.path.join(output_dir, subgroup), exist_ok=True)
            with gzip.open(input_filepath, 'rb') as f_in, open(decompressed_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

            # Lire la première ligne pour obtenir l'accession
            with open(decompressed_path, "r", encoding='utf-8') as reader:
                # first_line = reader.readline().strip()
                first_line = reader.readline()
                accession = first_line.split('.')[0][1:]  # Extraction de l'accession
                accession_map[accession] = decompressed_path
                decompressed_files.append(decompressed_path)

            if decompressed_path is None:
                print(f"Warning: Accession {accession} not found in accession_map {decompressed_files}")

        except Exception as e:
            print(f"Error decompressing {genome_file}: {e}")

    # Étape 2 : Récupération des informations de taxonomie par batch
    accessions = list(accession_map.keys())
    print(accessions)
    import numpy as np
    A = np.array(accessions)
    print(f"Nb input acessions {len(accessions)} , unique {len(np.unique(A))}")
    for i in tqdm(range(0, len(accessions), batch_size), desc="Fetching taxonomy data"):
        batch = accessions[i:i + batch_size]

        # try:
        #     with Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text") as taxo_handle:
        #         records = SeqIO.parse(taxo_handle, 'genbank')

        try:
            with Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text") as taxo_handle:
                # records = SeqIO.read(taxo_handle, 'genbank')
                records = SeqIO.parse(taxo_handle, 'genbank')
                # for record in records:
                #     accession = record.id.split('.')[0]
                #     classif = record.annotations.get('taxonomy', [])
                #     organism = record.annotations.get('organism', "Unknown Organism")
                for idx , record in enumerate(records):
                    accession = record.id.split('.')[0]  # Exclure la partie après le point
                    decompressed_path = accession_map.get(accession)
                    taxonomy = record.annotations.get('taxonomy', [])
                    organism = record.annotations.get('organism', "Unknown Organism")
                    print(f"Batch{i}, idx {idx }record.id {record.id} accession {accession} taxo len {len(taxonomy)}")
                    if not taxonomy:
                        print(f"Warning: No taxonomy found for accession {accession}.")
                        os.remove(decompressed_path)
                        continue

                    order = next((e for e in taxonomy if e.endswith('ales')), None)
                    group = taxonomy[2] if len(taxonomy) > 2 and not taxonomy[2].endswith('ales') else taxonomy[1]
                    if order:
                        # Construction du nom de fichier de sortie
                        file_name = f"{taxonomy[0]}_{taxonomy[1]}_{group}_{order}_{organism.replace(' ', '_')}.fna"
                        final_path = os.path.join(output_dir, subgroup, file_name)
                        print(f"order {order} group {group} ")
                        id_file = 0
                        temp_final_path = final_path
                        while os.path.exists(temp_final_path):
                            id_file +=1
                            temp_final_path =f"{final_path}_{id_file}"

                        final_path =temp_final_path
                        if id_file>0:
                            print(f"Target file rename with id_file {id_file} ")

                        try:
                            os.rename(decompressed_path, final_path)
                            print(f"Renamed: {decompressed_path} -> {final_path}")
                            nb_convert += 1
                        except Exception as e:
                            print(
                                f"Error renaming file {os.path.basename(decompressed_path)} -> {os.path.basename(final_path)}: {e}")
                    else:
                        print(f" no order {order} group {group} ")

                        # Suppression des fichiers sans ordre ales valide
                        os.remove(decompressed_path)
                        print(f"Removed {os.path.basename(decompressed_path)} (sans ordre ales)")


        except Exception as e:
                print(f"Error fetching or parsing data for accession {accession}: {e}")
                continue


    return  nb_convert
    # # Nettoyage des fichiers restants non traités
    # for file in decompressed_files:
    #     if os.path.exists(file):
    #         os.remove(file)
    #         print(f"Removed unprocessed file: {file}")

time_start = time.time()

datadir = "/home/hcourtei/Projects/MicroTaxo/codes/data/refseq"  # args.datadir "/groups/microtaxo/data/refseq"
subgroups = natsorted(os.listdir(datadir))
group_id = 1
batch_size = 10
nb_convert = get_taxo_by_batch(datadir, subgroups[group_id], batch_size)
duration = time.time() - time_start
print(f" Get Taxo for {datadir}, 'group_id' {group_id}")
print(f" nb_convert in output dir : {nb_convert}")