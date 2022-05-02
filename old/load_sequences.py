from python_tools import my_parser, my_pretty_printer
import os


def load_sample(filename: str, size_kmer) -> list[Sample]:
    "Loads from a file a list of samples"
    data = my_parser(filename, clean=True, merge=False)
    return [Sample(
            seq_dna=value,
            kmer_size=size_kmer,
            seq_id=key,
            seq_specie=filename.split('/')[-1].split('.')[0])
            for key, value in data.items()]


def load_list_sample(path: str, kmer: int) -> list[Sample]:
    "Loads from a folder a list of samples"
    my_samples: list[Sample] = []
    lst_sequences = listdir(path)
    for e in lst_sequences:
        mn = f"{e.split('_')[1]}_{e.split('_')[2].split('.')[0]}"
        parsed = my_parser(f"{path}{e}", clean=True, merge=True,
                           merge_name=mn)
        my_samples.append(
            Sample(
                seq_dna=parsed[mn],
                seq_id=mn,
                seq_specie=e,
                kmer_size=kmer))
    return my_samples


def species_list(input_dir: str) -> list:
    """
    Return the list of all species in a given dir, based upon their filenames
    """
    return [file.split('.')[0] for file in os.listdir(input_dir)]


def species_map(input_dir: str) -> dict:
    """
    Return a dict, association of a species name and a code for species
    """
    species_listed = list(set([file.split('_')[0]
                          for file in species_list(input_dir)]))
    return {species_listed[i]: i for i in range(len(species_listed))}


def sequences(input_dir: str):
    """
    Returns a dict of dicts, with all sequences inside, identified by the file name
    """
    full_dict: dict = {}
    files_list = species_list(input_dir)
    for i in range(len(files_list)):
        new_datas: dict = my_parser(
            f"{input_dir}{files_list[i]}.fna", True, True, files_list[i])
        full_dict = dict(full_dict, **new_datas)
    return full_dict


def sequence(input_file: str):
    """
    Returns a dict of dicts, with all sequences inside, identified by the file name
    """
    full_dict: dict = {}
    files_list = [input_file]
    for i in range(len(files_list)):
        new_datas: dict = my_parser(
            f"{input_file}", True, True, input_file.split('/')[-1][:-3])
        full_dict = dict(full_dict, **new_datas)
    return full_dict


if __name__ == "__main__":
    seqs: dict = sequences("/udd/sidubois/Stage/StreptoThermoGenomes/")
    my_pretty_printer(seqs, 10)
    print(len(seqs))
