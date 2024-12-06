"""Download genomes from NCBI."""
from os import system
from argparse import ArgumentParser, SUPPRESS
from rich.traceback import install
from Bio import SeqIO, Entrez
import time

def download_from_ncbi(summary_file: str, genomes_path: str, email: str, start: int = 0) -> None:
    """Download all genomes from a file, assumming its a standard NCBI summary file

    Args:
        summary_file (str): Path to a NCBI ftp file
        genomes_path (str): Path to folder to output genomes
        start (int, optional): A line number where to start from. Defaults to 0.
    """
    Entrez.email = email
    nb = 0
    all_duration = list()
    time_start = time.time()
    with open(summary_file, "r", encoding='utf-8') as summary_reader:
        # Skipping the two first lines
        next(summary_reader)
        next(summary_reader)

        for i, line in enumerate(summary_reader):
            if i > start and 'complete genome' in line.lower():  # Only keeping representative genomes

                split = line.split()

                try:
                    https_wget = [
                        elt for elt in split if elt.startswith('https')][0]
                    access = [split[0], https_wget]
                    # Downloading genome file
                    system(
                        f"wget -P {genomes_path} {access[1][8:]}/{access[1][8:].split('/')[-1]}_genomic.fna.gz && gzip -d {genomes_path}/{access[1][8:].split('/')[-1]}_genomic.fna.gz")
                    # Extracting
                    with open((file_path := f"{genomes_path}/{access[1][8:].split('/')[-1]}_genomic.fna"), "r", encoding='utf-8') as reader:
                        accession: str = reader.readline().split('.')[0][1:]

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
                            file_name: str = f"{genomes_path}/{classif[0]}_{classif[1]}_{group}_{order}_{sub.split(' ')[0]}_{sub.split(' ')[1]}.fna"
                            system(f"mv {file_path} {file_name}")
                        else:
                            # we clean genomes we can't retrive classification for
                            system(f"rm {file_path}")
                except Exception as exc:
                    raise BaseException(f"Job stopped at line {i}.") from exc


if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "summary_file", type=str, help="Path to NCBI summary ftp file.")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Download and rename a set of reference genomes')
    parser.add_argument(
        "-o", "--output", help="Specify a output folder", required=True)
    parser.add_argument(
        "-e", "--email", help="Specify a valid email for Entrez", required=True)
    parser.add_argument("-s", "--start", help="Specify a line number in your summary file to start from",
                        required=False, type=int, default=0)
    args = parser.parse_args()

    install(show_locals=True)

    # system(f"wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -o {args.output}")

    download_from_ncbi(args.summary_file, args.output,
                       args.email, start=args.start)
