"Functions to generate params file"

from json import load, dump
from argparse import ArgumentParser


def load_json(json_file: str) -> dict:
    """
    Charge un fichier json en un dictionnaire
    * json_file (str) : le chemin d'accès au fichier
    """
    path_to_json = f"parameters_files/{json_file.split('.')[0]}.json" if 'parameters_files' not in json_file else f"{json_file.split('.')[0]}.json"
    return load(open(path_to_json, "r"))


def save_json(json_file: str, dico_save: dict) -> None:
    """
    Sauve un dictionnaire en un fichier json
    * json_file (str) : le chemin d'accès au fichier
    * dico_save (dict) : le dictionnaire à sauvegarder
    """
    dump(dico_save, open(f"{json_file}.json", "w"), indent=4)


def my_params(filename: str):
    sparkle: str = "default"
    # dict to be converted in .json file to create our parameters set
    params_job: dict = {
        # subreads size and max. of reads in sample
        'window_size': 10000,
        'sampling_objective': 500,
        # params for your database here [kmer_size, subsampling_depth, pattern]
        'domain_ref': [5, 50, [1, 1, 1, 1, 1]],
        'phylum_ref': [5, 100, [1, 1, 1, 1, 1]],
        'group_ref': [4, 100, [1, 1, 1, 1]],
        'order_ref': [4, 100, [1, 1, 1, 1]],
        'family_ref': [4, 100, [1, 1, 1, 1]],
        'merged_ref': [4, 50, [1, 1, 1, 1]],
        # params for your sample here [kmer_size, pattern]
        'domain_sample': [5, [1, 1, 1, 1, 1]],
        'phylum_sample': [5, [1, 1, 1, 1, 1]],
        'group_sample': [4, [1, 1, 1, 1]],
        'order_sample': [4, [1, 1, 1, 1]],
        'family_sample': [4, [1, 1, 1, 1]],
        'merged_sample': [4, [1, 1, 1, 1]],
        # 'input' : location of reference genomes
        'input_train': f"/scratch/sdubois/lb_minion/train_{sparkle}/",
        # 'input_unk' : location of unk genomes
        'input_unk': f"/scratch/sdubois/lb_minion/unk_{sparkle}/",
        # 'output' : output for database
        'database_output': f"/scratch/sdubois/lb_minion/data/",
        'reports_output': f"/scratch/sdubois/lb_minion/output/{sparkle}/",
        # parameters for exploration and algorithm
        'threshold': 0.10,
        'nb_boosts': 10,
        'tree_depth': 10,
        # parameters regarding results
        'test_mode': 'no_test',  # 'no_test', 'min_set', 'verbose'
        # parameter for read selection, significance for softprob
        'reads_th': 0.1,
        'selection_mode': 'delta_mean',  # 'min_max','delta_mean','delta_sum'
        # force rebuilding full model, when you re-use a database but you changed model parameters
        'force_model_rebuild': False,  # never set true in multithread mode
        # tells the software if should consider both orientations or assume it is 3' -> 5' and computes canonical kmers
        'single_way': True,
        'targeted_level': 'family',  # domain, phylum, group, order, family
        'levels_list': ['domain', 'phylum', 'group', 'order', 'family'],
        'abundance_threshold': 0.25,
        # to fetch WISP genomes from refseq
        'email': 'siegfried.dubois@inria.fr',
        'annotate_path': 'genomes/to_annotate',
        'accession_numbers': 'genomes/assembly_summary.txt',
        # name for database
        'db_name': sparkle,
        'prefix_job': sparkle,
        'log_file': f'LOG_wisp_{sparkle}'
    }
    save_json(f"parameters_files/{filename}", params_job)


if __name__ == "__main__":
    parser = ArgumentParser()
    # declaring args
    parser.add_argument(
        "name", help="name for parameters file", type=str)
    # executing args
    args = parser.parse_args()
    my_params(args.name)
