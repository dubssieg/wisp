"""get taxonomy name for genomes from NCBI."""
import os
from Bio import SeqIO, Entrez
from argparse import ArgumentParser
# from multiprocessing import Pool
from tqdm import tqdm
import logging
from functools import partial
# import time
import gzip
import shutil




def get_taxo_name_refseq(input_dir: str , subgroup: str) -> None:
    """get taxonomy name for genomes from NCBI
    Args:
        genomes_dir (str): Path to folder with raw genomes xxx.fna.gz
    """
    output_dir = input_dir + '_unzip_with_taxo'

    Entrez.email = "hermann.courteille@irisa.fr"
    input_dir_sub= os.path.join(input_dir, subgroup)
    genom_file_list = os.listdir(input_dir_sub)

    for genom_file in tqdm(genom_file_list):

            input_filepath = os.path.join(input_dir_sub, genom_file)

            if input_filepath.endswith('.gz'):

                # Décompression
                try:
                    new_file_path = input_filepath[:-3] #drop .gz
                    output_file = os.path.join(output_dir, subgroup, os.path.basename(new_file_path))
                    os.makedirs(os.path.join(output_dir, subgroup), exist_ok=True)
                    return_code = os.system(f"gzip -dc {input_filepath} > {output_file} ")
                    if return_code != 0:
                        logger.error(f"Error during decompression of {os.path.basename(input_filepath)}")
                        continue  # Skip this file if decompression failed
                    # OU
                    with gzip.open(input_filepath, 'rb') as f_in:
                        with open(output_file, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    # # logger.info(f"Fichier décompressé avec succès : {output_file}")
                except Exception as e:
                    logger.error(f"Erreur lors de la décompression : {e}")

                # return_code = os.system(f"gzip -d {file_path} > /dev/null 2>&1")

                # if return_code == 0:
                #     # logger.info(f"Unzipped successfully: {file_path}")
                #     file_path = file_path[:-3]  # drop .gz
                #     time.sleep(0.2)
                # else:
                #     logger.error(f"Error during unzipping {os.path.basename(file_path)}. Return code: {return_code}")
                #     continue  # Skip file if unzip fails
            else:
                logger.warning(f"SKIP {os.path.basename(input_filepath)} not ending with .gz")
                continue

            try:

                with open(output_file, "r", encoding='utf-8') as reader:
                    line = reader.readline()
                    accession: str = line.split('.')[0][1:]

                # Get taxonomy from NCBI taxonomy
                with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as taxo_handle:
                    x = SeqIO.read(taxo_handle, 'genbank')
                    classif = x.annotations['taxonomy']
                    sub = x.annotations['organism']
                    order: int | str = 0
                    for e in classif:
                        if e[-4:] == 'ales':
                            order = e
                    if order:
                        group = classif[2] if classif[2][-4:] != 'ales' else classif[1]
                        file_name: str = f"{output_dir}/{subgroup}/{classif[0]}_{classif[1]}_{group}_{order}_{sub.split(' ')[0]}_{sub.split(' ')[1]}.fna"
                        os.system(f"mv {output_file} {file_name}")
                        # logger.info(f"Renamed {os.path.basename(file_path)} to {os.path.basename(file_name)}")
                    else:
                        # we clean genomes we can't retrive classification for
                        # os.system(f"rm {file_path}")
                        logger.warning(f" File {os.path.basename(input_filepath)}: classification not found")

            except Exception as e:
                logger.error(f"Error Accessing taxonomy {os.path.basename(input_filepath)}: {e}")
                # os.system(f"rm {file_path}")



if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    # parser.add_argument(
    #     "-d", "--datadir", help="Specify a output folder", required=True)
    parser.add_argument(
        "-i", "--id_group", help="Specify a output folder", required=True, type=int, default=0)
    args = parser.parse_args()
    id_group = args.id_group
    datadir = "/groups/microtaxo/data/refseq"
              #"/home/hcourtei/Projects/MicroTaxo/codes/data/refseq"  # args.datadir


    log_file = f"{datadir}/taxonomy_downloader{id_group}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)
    subgroups = sorted(os.listdir(datadir))
    logger.info(f"subdir {subgroups}")
    # max_processes = 4  # Nombre de processus parallèles
    # logger.info(f"Starting processing with {max_processes} processes for {len(subgroups)} subgroups")

    process_with_fixed_dir = partial(get_taxo_name_refseq, datadir)
    # for id_group, subgroup in enumerate(subgroups):

    logger.info('-' * 30)
    logger.info(f" Group {id_group}/{len(subgroups)} , Get Taxonomy from NCBI for in {datadir}/{subgroups[id_group]}")
    process_with_fixed_dir(subgroups[id_group])


    # with Pool(processes=max_processes) as pool:
    #     list(tqdm(pool.imap(process_with_fixed_dir, subgroups), total=len(subgroups)))

    logger.info("Processing completed.")
    # for id_group, subgroup in enumerate(subgroups):
    #     print('-'*30)
    #     # id_group, subgroup = 0, subgroups[0]
    #     print(f" Group {id_group}/{len(subgroups)} , Get Taxonomy from NCBI for in {datadir}")
    #     genome_dir = os.path.join(datadir, subgroup)
    #     get_taxo_name_refseq(genome_dir)
