"""get taxonomy name for genomes from NCBI."""
import os
from Bio import SeqIO, Entrez
from argparse import ArgumentParser
import tqdm


def get_taxo_name_refseq(genomes_dir: str) -> None:
    """get taxonomy name for genomes from NCBI
    Args:
        genomes_dir (str): Path to folder with raw genomes xxx.fna.gz
    """
    Entrez.email = "hermann.courteille@irisa.fr"
    genom_file_list = os.listdir(genomes_dir)


    for genom_file in tqdm.tqdm(genom_file_list):

            file_path = os.path.join(genomes_dir, genom_file)
            if file_path.endswith('.gz'):
                os.system(f"gzip -d {file_path}")
                file_path = file_path[:-3] # drop .gz
            else:
                print(f"SKIP {file_path} not ending with .gz")
                continue
            with open(file_path, "r", encoding='utf-8') as reader:
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
                    file_name: str = f"{genomes_dir}/{classif[0]}_{classif[1]}_{group}_{order}_{sub.split(' ')[0]}_{sub.split(' ')[1]}.fna"
                    os.system(f"mv {file_path} {file_name}")
                else:
                    # we clean genomes we can't retrive classification for
                    os.system(f"rm {file_path}")



if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "-d", "--datadir", help="Specify a output folder", required=True)

    args = parser.parse_args()

    get_taxo_name_refseq(args.datadir)
