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
        'domain_ref': [4, 10000, 30],
        'phylum_ref': [4, 10000, 400],
        'group_ref': [4, 10000, 300],
        'order_ref': [4, 10000, 300],
        'family_ref': [4, 10000, 200],
        # params for your sample here
        'domain_sample': [4, 10000, 200],
        'phylum_sample': [4, 10000, 600],
        'group_sample': [4, 10000, 600],
        'order_sample': [4, 10000, 600],
        'family_sample': [4, 10000, 500],
        # 'input' : location of genomes
        'input': "/udd/sidubois/Stage/Genomes/",
        # 'output' : output for database
        'output': "data/",
        # parameters for exploration and algorithm
        'threshold': 0.1,
        'nb_boosts': 10,
        'tree_depth': 8,
        # parameters regarding results
        'full_test_set': False,
        # parameter for read selection, signifiance for softprob
        'reads_th': 0.1,
        'selection_mode': 'delta_mean',  # 'min_max','delta_mean','delta_sum'
        # force rebuilding full model, when you re-use a database but you changed model parameters
        'force_model_rebuild': True
    }
    save_json("wisp_params", params_job)


my_params()
