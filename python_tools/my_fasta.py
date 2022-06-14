from Bio import SeqIO, Entrez
from os import path


def my_parser(filename: str, clean: bool = False, merge: bool = False, merge_name: str = "Merged") -> dict:
    """
    Renvoie un dictionnaire contenant toutes les séquences
    key : desc de la séquence
    value : séquence

    * filename le chemin de fichier
    * clean si le fichier doit être nettoyé de ses N
    * merge si on merge tous les fasta d'un fichier en une seule chaine
    """
    if clean and merge:
        loading = {fasta.id: str(fasta.seq).replace('N', '')
                   for fasta in SeqIO.parse(open(filename), 'fasta')}
        return {str(merge_name): ''.join([seq for seq in loading.values()])}
    elif clean:
        return {fasta.id: str(fasta.seq).replace('N', '') for fasta in SeqIO.parse(open(filename), 'fasta')}
    elif merge:
        loading = {fasta.id: str(fasta.seq)
                   for fasta in SeqIO.parse(open(filename), 'fasta')}
        return {str(merge_name): ''.join([seq for seq in loading.values()])}
    else:
        return {fasta.id: str(fasta.seq) for fasta in SeqIO.parse(open(filename), 'fasta')}


def my_fasta_parser(filename: str) -> dict:
    """Loads fasta files

    Args:
        filename (str): path to file

    Returns:
        dict: all sequences inside fasta file
    """
    return {fasta.id: str(fasta.seq).replace('N', '') for fasta in SeqIO.parse(open(filename), 'fasta')}


def my_pretty_printer(seq_dict: dict, size: int = 10) -> None:
    for key, value in seq_dict.items():
        print(
            f"{key} -> [{value[:size]} --- {value[-size:]}] : {len(value)} bp")


def my_classification_mapper(file):
    Entrez.email = ""  # EMAIL
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
        except:
            print(f"Can't get data for {file}")
            return None
    return None


def my_fetcher(filelist, outname):
    """
    Fetches all seq to one single file
    """
    Entrez.email = 'siegfried.dubois@inria.fr'
    for file in filelist:

        if not path.isfile(f"gen/{file}_{outname}.fna"):

            with Entrez.efetch(db="nucleotide", id=file,
                               rettype="fasta", retmode="text") as handle:
                with open(f"gen/{outname}.fna", "w") as writer:
                    writer.write(handle.read())
