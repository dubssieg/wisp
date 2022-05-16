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
        # 'taxa' : [kmer_size, reads_size, subsampling_depth]
        # params for your database here
        'domain_ref': [4, 10000, 10],
        'phylum_ref': [4, 10000, 150],
        'group_ref': [4, 10000, 150],
        'order_ref': [4, 10000, 150],
        'family_ref': [4, 10000, 100],
        # params for your sample here
        'domain_sample': [4, 10000, 200],
        'phylum_sample': [4, 10000, 1000],
        'group_sample': [4, 10000, 1000],
        'order_sample': [4, 10000, 1000],
        'family_sample': [4, 10000, 500],
        # 'input' : location of genomes
        'input': "/udd/sidubois/Stage/Genomes/",
        # 'output' : output for database
        'output': "data/",
        # parameters for exploration and algorithm
        'threshold': 0.2,
        'nb_boosts': 10,
        'tree_depth': 8,
        # parameters regarding results
        'full_test_set': False,
        # parameter for read selection, signifiance for softprob
        'reads_th': 0.25,
        'selection_mode': 'delta_mean',  # 'min_max','delta_mean','delta_sum'
        # force rebuilding full model, when you re-use a database but you changed model parameters
        'force_model_rebuild': True
    }
    save_json("wisp_params", params_job)


my_params()
