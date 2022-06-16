"This script is rather a set of tools for downloading references ganomes and renaming those"

from python_tools import my_classification_mapper, my_fetcher, my_parser, my_types_checker
from os import listdir, rename, system
import csv
from Bio import SeqIO
from random import random, choice
from argparse import ArgumentParser


def rename_genomes(path_to_genomes: str) -> None:
    """Aims to rename each genome it finds in the folder according to the WISP system.
    Only works on .fna files, as its the default input format for WISP.

    Args:
        path_to_genomes (str, optional): path to files to rename. Defaults to ANNOTATE_PATH.
    """
    files = [f"{file[:-4]}" for file in listdir(
        f"{path_to_genomes}/") if file[-4:] == '.fna']

    for file in files:
        new_name = my_classification_mapper(
            file, 'siegfried.dubois@inria.fr')
        print(new_name)

        if new_name != None:
            rename(f"{path_to_genomes}/{file}.fna",
                   f"{path_to_genomes}/{new_name}.fna")


def pre_rename(annotate_path: str) -> None:
    files = listdir(f"{annotate_path}/")

    for file in files:
        if('Bacteria' not in file and 'Archaea' not in file):
            with open(f"{annotate_path}/{file}", "r") as reader:
                accession = reader.readline().split('.')[0][1:]
            rename(f"{annotate_path}/{file}",
                   f"{annotate_path}/{accession}.fna")


def get_info(csv_file: str) -> None:
    with open(csv_file) as file:
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows = [[row[header.index("GenBank ID/Accession")], row[header.index(
            "Domain")], row[header.index("Phylum")], row[header.index("Class")], row[header.index("Order")], row[header.index("Family")], f"{row[header.index('Genome Name')].split(' ')[0]}_{row[header.index('Genome Name')].split(' ')[1]}"] for row in csvreader]

    for i, row in enumerate(rows):
        try:
            my_fetcher(row[0].split(',')[0], '_'.join(row[1:]))
        except:
            print("Bad request")


def gather(csv_file: str) -> None:
    ids = []
    with open(csv_file) as file:
        csvreader = csv.reader(file)
        header = next(csvreader)
        for row in csvreader:
            ids = list(
                ids + row[header.index("GenBank ID/Accession")].split(', '))
    with open("out.txt", "w") as writer:
        for id in ids:
            writer.write(f"{id}\n")

# downloading via https://www.ncbi.nlm.nih.gov/sites/batchentrez
# retrieving one huge .fasta file


def scatter(fasta_file: str) -> None:
    my_tuple = [(fasta.id, fasta.seq)
                for fasta in SeqIO.parse(open(fasta_file), 'fasta')]
    for id, seq in my_tuple:
        name: str | None = id  # my_classification_mapper(id)
        if name != None:
            with open(f"gen/{name}.fna", "w") as writer:
                writer.write(f"> {id} {name.replace('_',' ')}\n{seq}")
            print(f"File {name}.fna sucessfully writed out!")


def verificator(fasta_file: str) -> int:
    """Counts number of fasta sequences inside a huge fasta file

    Args:
        fasta_file (str): a path to a file

    Returns:
        int: number of subsequences
    """
    integer = 0
    with open(fasta_file, "r") as reader:
        for line in reader:
            if '>' in line:
                integer += 1
    return integer


def summary_to_dl(summary_file: str) -> None:
    # assumming its a standard NCBI summary file
    with open(summary_file, "r") as summary_reader:
        next(summary_reader)
        next(summary_reader)
        for i, line in enumerate(summary_reader):
            if i < 60000:
                next(summary_reader)
            elif i >= 60000:
                # accession, https
                split = line.split()
                print(f"Resolving entry n°{i} : {split[0]}")
                try:
                    https_wget = [
                        elt for elt in split if elt[:5] == 'https'][0]
                    access = [split[0], https_wget]
                    system(
                        f"wget -P {ANNOTATE_PATH} {access[1][8:]}/{access[1][8:].split('/')[-1]}_genomic.fna.gz; gzip -d {ANNOTATE_PATH}/*.gz")
                    pre_rename()
                    rename_genomes()
                except:
                    pass
            else:
                break


def destroy_sequence(sequence_path: str, destruction_ratio: float) -> None:
    for seq in listdir(sequence_path):
        header, sequence = "", ""
        with open(f"{sequence_path}{seq}", 'r') as reader:
            for line in reader:
                if line[0] == '>':
                    header = line
                else:
                    sequence = f"{sequence}{line}"
        new_seq = ''.join([base if random() >= destruction_ratio else choice(
            ['A', 'T', 'C', 'G']) for base in sequence])
        with open(f"{sequence_path}{seq.split('.')[0]}_destr.fna", 'w') as writer:
            writer.write(f"{header}\n{new_seq}")


def read_generator(path_to_sample: str, size_to_sample: int, start: int) -> None:
    "Sample one read out of each file in specified dir"
    size_to_sample, start = int(size_to_sample), int(start)
    files = listdir(path_to_sample)
    lectures = {file: my_parser(f"{path_to_sample}/{file}", True, True, file)[file]
                for file in files}
    reads = {k: v[start:start+size_to_sample]
             for k, v in lectures.items() if len(v) > size_to_sample+start}
    with open("sequences.fna", "w") as writer:
        writer.write('\n'.join(['> '+k+'\n'+v for k, v in reads.items()]))


def clean_rename(genomes_path: str) -> None:
    files = listdir(f"{genomes_path}/")

    for file in files:
        new_name = file
        for char in ['[', ']', '(', ')']:
            new_name = new_name.replace(char, '')
        rename(f"{genomes_path}/{file}",
               f"{genomes_path}/{new_name}")


if __name__ == "__main__":
    # This class aims to help to generate databases from a set of references.
    # From an assembly_summary.txt downloaded from NCBI (obtained via ftp.ncbi)
    # it downloads all references, renames it according to the WISP file system
    # All you have to do after is to put those files into train folder once process is done
    """
    parser = ArgumentParser()
    parser.add_argument(
        "-m", "--method", type=str, help="A callable func to execute")
    parser.add_argument('kwargs', nargs='*')
    args = parser.parse_args()
    try:
        getattr(tools, args.method)(*args.kwargs)
    except Exception as exc:
        raise BaseException("Bad execution") from exc
    """
    rename_genomes('genomes/143_prokaryote_genomes')
