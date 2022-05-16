from python_tools import my_classification_mapper, my_fetcher
from os import listdir, rename, system
import csv
from Bio import SeqIO
from constants import ANNOTATE_PATH, SUMMARY_FILE


def rename_genomes():
    files = [f"{file[:-4]}" for file in listdir(
        f"{ANNOTATE_PATH}/")]

    for file in files:
        new_name = my_classification_mapper(file)

        if new_name != None:
            rename(f"{ANNOTATE_PATH}/{file}.fna",
                   f"{ANNOTATE_PATH}/{new_name}.fna")


def cancel_genomes():
    files = [file for file in listdir(
        "/udd/sidubois/Stage/Genomes/train/") if len(file.split('_')) != 6]
    print(files)

    for file in files:
        new_name = f"_{file}"

        if new_name != None:
            rename(f"{ANNOTATE_PATH}/{file}.fna",
                   f"{ANNOTATE_PATH}/{new_name}.fna")


def pre_rename():
    files = listdir(f"{ANNOTATE_PATH}/")

    for file in files:
        if('Bacteria' not in file and 'Archaea' not in file):
            with open(f"{ANNOTATE_PATH}/{file}", "r") as reader:
                accession = reader.readline().split('.')[0][1:]
            rename(f"{ANNOTATE_PATH}/{file}",
                   f"{ANNOTATE_PATH}/{accession}.fna")


def get_info(csv_file):
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


def gather(csv_file):
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


def scatter(fasta_file):
    my_tuple = [(fasta.id, fasta.seq)
                for fasta in SeqIO.parse(open(fasta_file), 'fasta')]
    for id, seq in my_tuple:
        name: str | None = id  # my_classification_mapper(id)
        if name != None:
            with open(f"gen/{name}.fna", "w") as writer:
                writer.write(f"> {id} {name.replace('_',' ')}\n{seq}")
            print(f"File {name}.fna sucessfully writed out!")


def verificator(fasta_file):
    integer = 0
    with open(fasta_file, "r") as reader:
        for line in reader:
            if '>' in line:
                integer += 1
    print(integer)


def summary_to_dl(summary_file):
    # assumming its a standard NCBI summary file
    with open(summary_file, "r") as summary_reader:
        next(summary_reader)
        next(summary_reader)
        for i, line in enumerate(summary_reader):
            if i < 5500:
                print("Already done")
                next(summary_reader)
            elif i >= 5500 and i < 6500:
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


if __name__ == "__main__":
    # This class aims to help to generate databases from a set of references.
    # From an assembly_summary.txt downloaded from NCBI (obtained via ftp.ncbi)
    # it downloads all references, renames it according to the WISP file system
    # All you have to do after is to put those files into train folder once process is done
    summary_to_dl(SUMMARY_FILE)
