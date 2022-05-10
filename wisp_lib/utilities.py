# to do random tests

from python_tools import my_classification_mapper, my_fetcher
from wisp_lib import species_list
from os import listdir, rename
import csv
from Bio import SeqIO


def rename_genomes():
    files = [f"{file[:-4]}" for file in listdir(
        "/udd/sidubois/Stage/Genomes/143_genomes/")]

    for file in files:
        new_name = my_classification_mapper(file)

        if new_name != None:
            rename(f"/udd/sidubois/Stage/Genomes/143_genomes/{file}.fna",
                   f"/udd/sidubois/Stage/Genomes/143_genomes/{new_name}.fna")


def cancel_genomes():
    files = [file for file in listdir(
        "/udd/sidubois/Stage/Genomes/train/") if len(file.split('_')) != 6]
    print(files)

    for file in files:
        new_name = f"_{file}"

        if new_name != None:
            rename(f"/udd/sidubois/Stage/Genomes/143_genomes/{file}.fna",
                   f"/udd/sidubois/Stage/Genomes/143_genomes/{new_name}.fna")


def pre_rename():
    files = listdir("/udd/sidubois/Stage/Genomes/143_genomes/")

    for file in files:
        with open(f"/udd/sidubois/Stage/Genomes/143_genomes/{file}", "r") as reader:
            accession = reader.readline().split('.')[0][1:]
        rename(f"/udd/sidubois/Stage/Genomes/143_genomes/{file}",
               f"/udd/sidubois/Stage/Genomes/143_genomes/{accession}.fna")


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
        name: str | None = my_classification_mapper(id)
        if name != None:
            with open(f"gen/{name}.fna", "w") as writer:
                writer.write(f"> {id} {name.replace('_',' ')}\n{seq}")
            print(f"File {name}.fna sucessfully writed out!")


cancel_genomes()
