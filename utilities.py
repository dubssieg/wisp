"This script is rather a set of tools for downloading references genomes and renaming those"

from os import listdir, rename, system
import csv
from random import random, choice
from argparse import ArgumentParser
from Bio import SeqIO
from python_tools import my_classification_mapper, my_parser, my_output_msg, my_logs_global_config


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
        if new_name is not None:
            rename(f"{path_to_genomes}/{file}.fna",
                   f"{path_to_genomes}/{new_name}.fna")


def pre_rename(annotate_path: str) -> None:
    """Gets the accession number to later retrieve taxonomy, and puts it as filename

    Args:
        annotate_path (str): path we want to process
    """
    files = listdir(f"{annotate_path}/")

    for file in files:
        if('Bacteria' not in file and 'Archaea' not in file):
            with open(f"{annotate_path}/{file}", "r") as reader:
                accession = reader.readline().split('.')[0][1:]
            rename(f"{annotate_path}/{file}",
                   f"{annotate_path}/{accession}.fna")


def gather(csv_file: str) -> None:
    """Extacts accession numbers from CSV file

    Args:
        csv_file (str): a path to a .csv file
    """
    ids = []
    with open(csv_file) as file:
        csvreader = csv.reader(file)
        header = next(csvreader)
        for row in csvreader:
            ids = list(
                ids + row[header.index("GenBank ID/Accession")].split(', '))
    with open("out.txt", "w") as writer:
        for id_file in ids:
            writer.write(f"{id_file}\n")


def scatter(fasta_file: str) -> None:
    """When downloading via https://www.ncbi.nlm.nih.gov/sites/batchentrez
    we are retrieving one huge .fasta file ; this func splits it in independant files

    Args:
        fasta_file (str): a huge .fna file
    """
    my_tuple = [(fasta.id, fasta.seq)
                for fasta in SeqIO.parse(open(fasta_file), 'fasta')]
    for id_seq, seq in my_tuple:
        name: str | None = id_seq
        if name is not None:
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
    """Download all genomes from a file, assumming its a standard NCBI summary file

    Args:
        summary_file (str): path to a NCBI file
    """
    genomes_path = "genomes/to_annotate"
    with open(summary_file, "r") as summary_reader:
        next(summary_reader)
        next(summary_reader)
        for i, line in enumerate(summary_reader):
            if i < 165000:
                next(summary_reader)
            elif i >= 165000:
                # accession, https
                split = line.split()
                my_output_msg(f"Resolving entry n°{i} : {split[0]}")
                try:
                    https_wget = [
                        elt for elt in split if elt[:5] == 'https'][0]
                    access = [split[0], https_wget]
                    system(
                        f"wget -P {genomes_path} {access[1][8:]}/{access[1][8:].split('/')[-1]}_genomic.fna.gz; gzip -d {genomes_path}/*.gz")
                    pre_rename(genomes_path)
                    rename_genomes(genomes_path)
                    # we clean genomes we can't retrive classification for
                    system(f"rm {genomes_path}/N*")
                    clean_rename(genomes_path)
                except OSError:
                    pass
            else:
                break


def destroy_sequence(sequence_path: str, sequence_output: str, destruction_ratio: float) -> None:
    """Applies artificial sequencing errors to sequences

    Args:
        sequence_path (str): path to reference sequences to be destroyed
        sequence_output (str): path where files will be outputted
        destruction_ratio (float): percentage we want to destroy our sequences
    """
    destruction_ratio = float(destruction_ratio)
    for i, seq in enumerate(listdir(sequence_path)):
        my_output_msg(f"Destroying sequence {i}...")
        header, sequence = "", ""
        with open(f"{sequence_path}{seq}", 'r') as reader:
            for line in reader:
                if line[0] == '>':
                    header = line
                else:
                    sequence = f"{sequence}{line}"
        new_seq = ''.join([base if random() >= 0.75*destruction_ratio else choice(
            ['A', 'T', 'C', 'G']) for base in sequence])
        with open(f"{sequence_output}{seq.split('.')[0]}_destr.fna", 'w') as writer:
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
    """Cleanses genomes from non-standard characters that might cause issues

    Args:
        genomes_path (str): path to apply cleaning
    """
    files = listdir(f"{genomes_path}/")

    for file in files:
        new_name = file
        for char in ['[', ']', '(', ')']:
            new_name = new_name.replace(char, '')
        for char in [' ']:
            new_name = new_name.replace(char, '-')
        rename(f"{genomes_path}/{file}",
               f"{genomes_path}/{new_name}")


if __name__ == "__main__":
    # This class aims to help to generate databases from a set of references.
    # From an assembly_summary.txt downloaded from NCBI (obtained via ftp.ncbi)
    # it downloads all references, renames it according to the WISP file system
    # All you have to do after is to put those files into train folder once process is done

    parser = ArgumentParser()
    parser.add_argument(
        "-m", "--method", type=str, choices=['clean_rename', 'summary_to_dl', 'destroy_sequence'], help="A callable func to execute")
    parser.add_argument('kwargs', nargs='*')
    args = parser.parse_args()

    my_logs_global_config("WISP_utilities", True, True)
    match args.method:
        case 'clean_rename':
            # Need to provide path to target folder
            clean_rename(**args.kwargs)
        case 'summary_to_dl':
            # Need to provide mpath to accession file
            summary_to_dl(**args.kwargs)
        case 'destroy_sequence':
            # input folder, output folder, destruction ratio
            destroy_sequence(**args.kwargs)
        case _:
            my_output_msg(
                'You need to specifiy a method and its args for the program to work. See documentation for help')
