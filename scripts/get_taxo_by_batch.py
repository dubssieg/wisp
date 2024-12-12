import os
import gzip
import shutil
import time
from natsort import natsorted
from Bio import Entrez, SeqIO
import logging
import numpy as np
from tqdm import tqdm
from argparse import ArgumentParser
import yaml

parser = ArgumentParser(add_help=False)
parser.add_argument( "-datadir", type=str, help="Path to refseq zipped data .gz")
parser.add_argument("-group_id", type=int,  help="Specify a output folder")
parser.add_argument("-batchsize", type=int,  help="Specify a output folder")

args = parser.parse_args()
# > python get_taxo_by_batch.py -datadir /groups/microtaxo/data/refseq -group_id 0 -batchsize 100


log_file = f"unzip_and_get_taxo_group_{args.group_id}.log"
if os.path.exists(log_file):
    os.remove(log_file)

logging.basicConfig(
    level=logging.INFO,  # Change to DEBUG for more detailed logs
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%H:%M",  # Format de l'heure : heures et minutes
    handlers=[logging.FileHandler(log_file),  # Logs to a file
              logging.StreamHandler()  # Logs to consol
             ]
)
logger = logging.getLogger(__name__)



def get_taxo_by_batch(input_dir: str, subgroup: str, batch_size: int = 100, api_key=None, mail=None, verbose=False) -> int:
    """
    Get taxonomy name for genomes from NCBI using batch requests.

    Args:
        input_dir (str): Path to folder with raw genome files (.fna.gz).
        subgroup (str): Subgroup name for processing specific folder.
        batch_size (int): Number of accessions to fetch in a single request.
    """
    # Configurer les paramètres Bio.Entrez
    Entrez.email = mail  # Remplacez par votre email
    Entrez.api_key = api_key #""
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
                logger.warning(f"Warning: Accession {accession} not found in accession_map {decompressed_files}")

        except Exception as e:
            logger.error(f"Error decompressing {genome_file}: {e}")

    # Étape 2 : Récupération des informations de taxonomie par batch
    accessions = list(accession_map.keys())
    if verbose:
        logger.info(f"Nb input acessions {len(accessions)} , unique {len(np.unique(np.array(accessions)))}")
    for i in tqdm(range(0, len(accessions), batch_size), desc="Fetching taxonomy data"):
        batch = accessions[i:i + batch_size]
        try:
            with Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text") as taxo_handle:
                # records = SeqIO.read(taxo_handle, 'genbank')
                records = SeqIO.parse(taxo_handle, 'genbank')
                for idx , record in enumerate(records):
                    accession = record.id.split('.')[0]  # Exclure la partie après le point
                    decompressed_path = accession_map.get(accession)
                    taxonomy = record.annotations.get('taxonomy', [])
                    organism = record.annotations.get('organism', "Unknown Organism")
                    if verbose:
                        logger.info(f"Batch {i}, idx {idx} record.id {record.id} taxo len {len(taxonomy)}")

                    if not taxonomy:
                        logger.warning(f"Warning: No taxonomy found for accession {accession}.")
                        os.remove(decompressed_path)
                        continue

                    order = next((e for e in taxonomy if e.endswith('ales')), None)

                    if order:
                        regne, phylum = taxonomy[0], taxonomy[1]
                        group = taxonomy[2] if len(taxonomy) > 2 and not taxonomy[2].endswith('ales') else taxonomy[1]

                        # Construction du nom de fichier de sortie
                        file_name = f"{regne}_{phylum}_{group}_{order}_{organism.replace(' ', '_')}.fna"
                        final_path = os.path.join(output_dir, subgroup, file_name)
                        # logger.info(f"order {order} group {group} ")
                        id_file = 1
                        temp_final_path = final_path
                        while os.path.exists(temp_final_path):
                            temp_final_path = f"{final_path.rsplit('.fna', 1)[0]}_{id_file}.fna"
                            id_file +=1

                        final_path = temp_final_path

                        try:
                            os.rename(decompressed_path, final_path)
                            if verbose:
                                logger.info(f"Renamed: {os.path.basename(decompressed_path)} -> {os.path.basename(final_path)}")
                            nb_convert += 1
                        except Exception as e:
                            logger.error(f"Error renaming file {os.path.basename(decompressed_path)} "
                                  f"-> {os.path.basename(final_path)}: {e}")
                    else:
                        logger.warning(f" no order {order} group {group} e.endswith('ales')")
                        # Suppression des fichiers sans ordre ales valide
                        os.remove(decompressed_path)
                        logger.warning(f"Removed {os.path.basename(decompressed_path)}")

        except Exception as e:
                logger.error(f"Error fetching or parsing data for accession {accession}: {e}")
                continue

    return  nb_convert

time_start = time.time()

datadir = args.datadir #"/home/hcourtei/Projects/MicroTaxo/codes/data/refseq"  # args.datadir "/groups/microtaxo/data/refseq"
subgroups = natsorted(os.listdir(datadir))
group_id = args.group_id #1
batchsize = args.batchsize
# Chemin vers le fichier YAML
yaml_file = "NCBI_credentials.yaml"

with open(yaml_file, 'r') as file:
    credentials = yaml.safe_load(file)

mail = credentials.get("mail")
api_key = credentials.get("api_key")

nb_convert = get_taxo_by_batch(datadir, subgroups[group_id], batchsize, api_key, mail, verbose=False)
duration = time.time() - time_start
logger.info(f" Get Taxo for {datadir}, 'group_id' {group_id}")
logger.info(f" {nb_convert} files with taxo, extract in {round(duration/60)} min")
# 19:13 - INFO -  5000 files with taxo, extract in 91 min
