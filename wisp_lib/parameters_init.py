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
        'domain_ref': [4, 10000, 50],
        'phylum_ref': [4, 10000, 250],
        'group_ref': [4, 10000, 375],
        'order_ref': [4, 10000, 375],
        'family_ref': [4, 10000, 250],
        # params for your sample here
        'domain_sample': [4, 10000, 500],
        'phylum_sample': [4, 10000, 750],
        'group_sample': [4, 10000, 750],
        'order_sample': [4, 10000, 750],
        'family_sample': [4, 10000, 500],
        # 'input' : location of genomes
        'input': "/udd/sidubois/Stage/Genomes/",
        # 'output' : output for database
        'output': "data/",
        # parameters for exploration and algorithm
        'threshold': 0.1,
        'nb_boosts': 10,
        # parameters regarding restuls
        'full_test_set': False
    }
    save_json("wisp_params", params_job)


my_params()
