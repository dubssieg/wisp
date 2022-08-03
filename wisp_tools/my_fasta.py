from os import path
from Bio import SeqIO, Entrez


def my_parser(filename: str, clean: bool = False, merge: bool = False, merge_name: str = "Merged", objective: int = 10000) -> dict[str, str]:
    """
    Renvoie un dictionnaire contenant toutes les séquences
    key : desc de la séquence
    value : séquence

    * filename le chemin de fichier
    * clean si le fichier doit être nettoyé de ses N
    * merge si on merge tous les fasta d'un fichier en une seule chaine
    """
    match clean, merge:
        case True, True:
            loading = {fasta.id: ''.join([letter for letter in str(fasta.seq) if letter in ['A', 'T', 'C', 'G']])
                       for fasta in SeqIO.parse(open(filename), 'fasta')}
            ret = {str(merge_name): ''.join([seq for seq in loading.values()])}
        case True, False:
            ret = {fasta.id: ''.join([letter for letter in str(fasta.seq) if letter in [
                                     'A', 'T', 'C', 'G']]) for fasta in SeqIO.parse(open(filename), 'fasta')}
        case False, True:
            loading = {fasta.id: str(fasta.seq)
                       for fasta in SeqIO.parse(open(filename), 'fasta')}
            ret = {str(merge_name): ''.join([seq for seq in loading.values()])}
        case _:
            ret = {fasta.id: str(fasta.seq)
                   for fasta in SeqIO.parse(open(filename), 'fasta')}
    return {k: v for k, v in ret.items() if len(v) >= objective}


def my_fasta_parser(filename: str, length: int) -> dict[str, str]:
    """Loads fasta files

    Args:
        filename (str): path to file

    Returns:
        dict: all sequences inside fasta file
    """
    reads = {fasta.id: ''.join([letter for letter in str(fasta.seq) if letter in [
                               'A', 'T', 'C', 'G']]) for fasta in SeqIO.parse(open(filename), 'fasta')}
    return {k: v for k, v in reads.items() if len(v) >= length}


def my_pretty_printer(seq_dict: dict, size: int = 10) -> None:
    """Serves to show if sequences are correctly loaded

    Args:
        seq_dict (dict): sequences
        size (int, optional): size of window to show. Defaults to 10.
    """
    for key, value in seq_dict.items():
        print(
            f"{key} -> [{value[:size]} --- {value[-size:]}] : {len(value)} bp")


def my_classification_mapper(file: str, email: str):
    """Fetches the classification from NCBI taxonomy using accession number

    Args:
        file (str): name of the file
        email (str): identity

    Returns:
        str|none: taxonomy or none
    """
    Entrez.email = email
    # skipping unnecessary calls for already processed files
    if('Bacteria' not in file and 'Archaea' not in file):
        try:
            handle = Entrez.efetch(db="nucleotide", id=file,
                                   rettype="gb", retmode="text")
            x = SeqIO.read(handle, 'genbank')  # def : genbank
            classif = x.annotations['taxonomy']
            sub = x.annotations['organism']
            has_order = False
            order = ""
            for e in classif:
                if e[-4:] == 'ales':
                    order = e
                    has_order = True
            if has_order:
                group = classif[2] if classif[2][-4:] != 'ales' else classif[1]
                return f"{classif[0]}_{classif[1]}_{group}_{order}_{sub.split(' ')[0]}_{sub.split(' ')[1]}"
            else:
                return None
        except Exception:
            return None


def my_minion(filename: str, output_minion: str):
    """Purges MinION from unecessary stuff

    Args:
        filename (str): input file
        output_minion (str): output file
    """
    header, sequence, collection = '', '', {}
    with open(filename, 'r', encoding="utf-8") as minion_reader:
        for i, line in enumerate(minion_reader):
            if i % 4 == 0:
                header = line.split(' ')[0]
            elif i % 4 == 1:
                sequence = line[:-1]
            else:
                collection[header] = sequence
    with open(f"{output_minion}/{(filename.split('/')[-1]).split('.')[0]}.fna", 'w', encoding="utf-8") as minion_writer:
        minion_writer.write(
            '\n'.join([f"> {k}\n{v}" for k, v in collection.items()]))


def my_fetcher(filelist: list[str], outname: str, email: str):
    """
    Fetches all seq to one single file
    """
    Entrez.email = email
    for file in filelist:

        if not path.isfile(f"gen/{file}_{outname}.fna"):

            with Entrez.efetch(db="nucleotide", id=file,
                               rettype="fasta", retmode="text") as handle:
                with open(f"gen/{outname}.fna", "w", encoding="utf-8") as writer:
                    writer.write(handle.read())
