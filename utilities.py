"This script is rather a set of tools for downloading references genomes and renaming those"

from os import listdir, rename, system
from csv import reader
from random import random, choice, randrange
from argparse import ArgumentParser, Namespace
from typing import Callable
from Bio import SeqIO
from wisp_tools import my_classification_mapper, my_parser, my_output_msg, my_logs_global_config, my_minion
from wisp_view import number_of_classes, compare, compdiff_plotting, plot_database_features
import matplotlib.pyplot as plt
from collections import Counter


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
        csvreader = reader(file)
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


def format_tool(genome_path):
    pre_rename(genome_path)
    rename_genomes(genome_path)
    clean_rename(genome_path)
    system(f"rm {genome_path}/N*")


def summary_to_dl(summary_file: str, genomes_path: str, start: int = 0) -> None:
    """Download all genomes from a file, assumming its a standard NCBI summary file

    Args:
        summary_file (str): path to a NCBI file
    """
    start = int(start)

    my_output_msg(f"Summary filepath : {summary_file}")

    with open(summary_file, "r") as summary_reader:
        next(summary_reader)
        next(summary_reader)
        for i, line in enumerate(summary_reader):
            if i > start:
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
                except Exception:
                    pass
                # we clean genomes we can't retrive classification for
                system(f"rm {genomes_path}/N*")
                clean_rename(genomes_path)
            else:
                my_output_msg(f"Skipping item {i}")


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
        for char in ['..']:
            new_name = new_name.replace(char, '.')
        rename(f"{genomes_path}/{file}",
               f"{genomes_path}/{new_name}")


def minion(genomes_path: str, output_minion: str) -> None:
    files = listdir(f"{genomes_path}/")
    for file in files:
        my_minion(f"{genomes_path}/{file}", output_minion)


def retrieve(tdir: str, output_folder: str):
    glbl = Counter()
    list_of_rf = [f"{tdir}/{d}" for d in listdir(tdir)]
    list_of_files = []
    for rf in list_of_rf:
        list_of_files.extend(
            [f"{rf}/{txtfile}" for txtfile in listdir(rf) if '.txt' in txtfile])
    for fil in list_of_files:
        with open(fil, 'r') as reader:
            glbl += Counter([line.split(':')[1].replace('\n', '')
                            for line in reader])
    with open(f"{output_folder}/output_{tdir.split('/')[-1]}.txt", 'w') as writer:
        writer.write('\n'.join([key.split('\'')[9]+":"+str(value)
                     for key, value in glbl.items()]))


def executor(func: Callable, argsm: list, unpack: bool, hstring: str) -> None:
    try:
        func(*argsm) if unpack else func(args)
    except:
        my_output_msg(hstring)


def create_mock_dataset(inpult_folder: str, output_folder: str, destruction_ratio=0) -> None:
    sequences = {}
    list_of_genomes = [
        f"{inpult_folder}/{file}" for file in listdir(inpult_folder)]
    for genome in list_of_genomes:
        parsed_genome = my_parser(genome, True, False, objective=300000)
        if parsed_genome != {}:
            range_select = randrange(150000, 300000)
            partial_sequence = (list(parsed_genome.values())[0])[
                0:range_select]
            sequences[f"> {genome.split('/')[-1][:-4]} (length:{range_select})"] = ''.join([base if random() >= 0.75*destruction_ratio else choice(
                ['A', 'T', 'C', 'G']) for base in partial_sequence])
    with open(f"{output_folder}/mock_dataset.fna", "w") as fna_writer:
        fna_writer.write(
            "\n".join([f"{k}\n{v}" for k, v in sequences.items()]))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "method", type=str, choices=['mock_dataset', 'format_tool', 'aggregate', 'database_features', 'kmers_signatures', 'compare_outputs', 'clean_rename', 'summary_to_dl', 'destroy_sequence', 'clean_minion', 'extract_genomes'], help="A callable func to execute")
    parser.add_argument('kwargs', nargs='*',
                        help="Args for Callable, see documentation for usage")
    args = parser.parse_args()

    my_logs_global_config("WISP_utilities", '0', True, True)
    plt.rcParams.update({'figure.max_open_warning': 0})

    match args.method:
        case 'mock_dataset':
            func, unpack, hstring = create_mock_dataset, True, "Func needs a genome directory name and a output directory name"
        case 'format_tool':
            func, unpack, hstring = format_tool, True, "Func needs a raw NCBI-downloaded genome folder"
        case 'aggregate':
            func, unpack, hstring = retrieve, True, "Func needs a WISP output directory name (relative or absolute path) and a output path for file"
        case 'database_features':
            func, unpack, hstring = plot_database_features, True, "Func needs a WISP database directory name (relative or absolute path)"
        case 'kmers_signatures':
            func, unpack, hstring = compdiff_plotting, True, "Func needs a directory name containing genomes (relative or absolute path)"
        case 'compare_outputs':
            func, unpack, hstring = compare, False, "Func needs a list of WISP output directories names (relative or absolute paths)"
        case 'clean_rename':
            func, unpack, hstring = clean_rename, True, "Func needs a directory name containing genomes (relative or absolute path)"
        case 'summary_to_dl':
            func, unpack, hstring = summary_to_dl, True, "Func needs path to a accession file, a output folder, and optionally a line number to strart from"
        case 'destroy_sequence':
            func, unpack, hstring = destroy_sequence, True, "Func needs a genome directory input path, and a empty directory output path (both relative or absolute paths) and a destruction ratio"
        case 'clean_minion':
            func, unpack, hstring = minion, True, "Func needs a directory name containing genomes (relative or absolute path)"
        case 'extract_genomes':
            func, unpack, hstring = number_of_classes, True, "Func needs a number of relatives, a genome directory input path, and a empty directory output path (both relative or absolute paths)"
        case _:
            func, unpack, hstring = exit, False, "No correct method was used. Exiting."
    executor(func, args.kwargs, unpack, hstring)
