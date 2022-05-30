"Functions to generate params file"

from json import load, dump


def load_json(json_file: str) -> dict:
    """
    Charge un fichier json en un dictionnaire
    * json_file (str) : le chemin d'accès au fichier
    """
    return load(open(f"{json_file.split('.')[0]}.json", "r"))


def save_json(json_file: str, dico_save: dict) -> None:
    """
    Sauve un dictionnaire en un fichier json
    * json_file (str) : le chemin d'accès au fichier
    * dico_save (dict) : le dictionnaire à sauvegarder
    """
    dump(dico_save, open(f"{json_file}.json", "w"), indent=4)


def my_params():
    # dict to be converted in .json file to create our parameters set
    params_job: dict = {
        # 'taxa' : [kmer_size, reads_size, subsampling_depth, pattern]
        # merge mode parameters
        'merged_ref': [4, 10000, 50, '1111'],
        'merged_sample': [4, 10000, 200, '1111'],
        # params for your database here
        'domain_ref': [4, 10000, 30, '1111'],
        'phylum_ref': [4, 10000, 600, '1111'],
        'group_ref': [4, 10000, 600, '1111'],
        'order_ref': [4, 10000, 400, '1111'],
        'family_ref': [4, 10000, 200, '1111'],
        # params for your sample here
        'domain_sample': [4, 10000, 200, '1111'],
        'phylum_sample': [4, 10000, 600, '1111'],
        'group_sample': [4, 10000, 600, '1111'],
        'order_sample': [4, 10000, 600, '1111'],
        'family_sample': [4, 10000, 500, '1111'],
        # 'input' : location of genomes
        'input': "/udd/sidubois/Stage/143_Genomes/",
        # 'output' : output for database
        'output': "data/",
        # parameters for exploration and algorithm
        'threshold': 0.10,
        'nb_boosts': 10,
        'tree_depth': 10,
        # parameters regarding results
        'full_test_set': True,
        # parameter for read selection, significance for softprob
        'reads_th': 0.10,
        'selection_mode': 'delta_mean',  # 'min_max','delta_mean','delta_sum'
        # force rebuilding full model, when you re-use a database but you changed model parameters
        'force_model_rebuild': False,  # never set true in multithread mode
        # tells the software if should consider both orientations or assume it is 3' -> 5' and computes canonical kmers
        'single_way': True,
        'targeted_level': 'family',  # domain, phylum, group, order, family
        'levels_list': ['domain', 'phylum', 'group', 'order', 'family']
    }
    save_json("wisp_params", params_job)


my_params()
