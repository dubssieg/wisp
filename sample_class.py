"This class is about the computing of DNA samples"
from typing import Callable
from wisp_tools import my_parser, my_function_timer
from wisp_lib import species_map, encoder, write_xgboost_data, load_mapping, splitting_generator, counter_ultrafast
from os import listdir
from json import dump
from pathlib import Path


def mapping_sp(input_dir: str, path: str, classif_level: str, db_name: str, int_level: int, taxa: str | None) -> dict:
    """Calls for a map of species, saves it, and return map as a dict

    Args:
        input_dir (str): dir where references genomes are stored
        path (str): path for database, used to save the map
        classif_level (str): _description_
        db_name (str): _description_
        int_level (int): _description_
        taxa (str | None): _description_

    Returns:
        dict: a map of species for gien taxa
    """
    list_sp = species_map(input_dir, int_level, taxa)
    my_path = f"{path}{db_name}/{classif_level}/{taxa}_saved_mapping.json" if taxa != None else f"{path}{db_name}/{classif_level}/saved_mapping.json"
    with open(my_path, 'w') as file_manager:
        dump(list_sp, file_manager)
    return list_sp


def mapping_merged_sp(input_dir: str, path: str, db_name: str) -> dict:
    """Calls for a map of species, saves it, and return map as a dict

    Args:
        input_dir (str): dir where references genomes are stored
        path (str): path for database, used to save the map
        classif_level (str): _description_
        db_name (str): _description_
        int_level (int): _description_
        taxa (str | None): _description_

    Returns:
        dict: a map of species for given taxa
    """
    list_sp = species_map(input_dir, 4)
    my_path = f"{path}{db_name}/merged/merged_saved_mapping.json"
    with open(my_path, 'w') as file_manager:
        dump(list_sp, file_manager)
    return list_sp


@my_function_timer("Building datasets")
def make_datasets(job_name: str, input_dir: str, path: str, datas: list[str], sampling: int, db_name: str, classif_level: str, kmer_size: int, read_size: int, pattern: list,  sp_determied: str | None):
    """
    Create the datasets and calls for storage

    * input_dir (str) : path to directory with reference genomes
    * path (str) : location where .txt.train, .txt.unk and .txt.test are stored
    * job_name (str) : name of job, identificator
    * classif_level (str) : level of cassification we're working at
    * db_name (str) : database we need to search in
    """
    my_encoder: dict = encoder(ksize=kmer_size)
    taxa: dict[str, int] = {
        'domain': 0,
        'phylum': 1,
        'group': 2,
        'order': 3,
        'family': 4
    }
    if datas == ['train', 'test']:
        # safe creation of dir
        Path(f"{path}{db_name}/{classif_level}/").mkdir(parents=True, exist_ok=True)
        # we re-generate database, so we need to map it out
        if classif_level != 'merged':
            sp_map = mapping_sp(input_dir, path,
                                classif_level, db_name, taxa[classif_level], sp_determied)
        else:
            sp_map = mapping_merged_sp(input_dir, path, db_name)
    else:
        # database already generated, loading mapping without erasing it
        sp_map = load_mapping(path, db_name,
                              classif_level, sp_determied)

    for type_data in datas:
        lst_sequences = [l for l in listdir(input_dir)] if classif_level == 'domain' or classif_level == 'merged' else [
            l for l in listdir(input_dir) if l.split('_')[taxa[classif_level]-1] == sp_determied]

        lines = []
        for e in lst_sequences:
            # sample
            genomes_scaffolds = my_parser(f"{input_dir}{e}", clean=True)

            for _, sequence in genomes_scaffolds.items():
                try:
                    # we create optimal list of reads
                    all_reads = splitting_generator(
                        sequence, read_size, int(sampling/3) if type_data == 'test' else sampling)
                    counters = [counter_ultrafast(split, kmer_size, pattern)
                                for split in all_reads]
                    encoded = [{my_encoder[k]:v for k, v in cts.items()}
                               for cts in counters]
                    lines = lines + \
                        [f"{sp_map[e.split('_')[taxa[classif_level]]]} {' '.join([str(k)+':'+str(v) for k,v in encoded_read.items()])} #{e}" for encoded_read in encoded]
                except Exception as exc:
                    raise exc

        if lines != []:
            write_xgboost_data(lines, path,
                               classif_level, type_data, db_name, job_name, sp_determied)


@my_function_timer("Building datasets")
def make_unk_datasets(all_reads, job_name: str, path: str, db_name: str, classif_level: str, kmer_size: int, pattern: list,  sp_determied: str | None, func: Callable) -> int:
    """
    Create the datasets and calls for storage

    * input_dir (str) : path to directory with reference genomes
    * path (str) : location where .txt.train, .txt.unk and .txt.test are stored
    * job_name (str) : name of job, identificator
    * classif_level (str) : level of cassification we're working at
    * db_name (str) : database we need to search in

    Returns number of subreads in dataset
    """
    my_encoder = encoder(kmer_size)
    counters = [counter_ultrafast(split, kmer_size, pattern)
                for split in all_reads]
    encoded = [{my_encoder[k]:v for k, v in cts.items()} for cts in counters]
    lines = [
        f"0 {' '.join([str(k)+':'+str(v) for k, v in counts.items()])}" for counts in encoded]
    write_xgboost_data(lines, path,
                       classif_level, 'unk', db_name, job_name, sp_determied)
    return len(lines)
