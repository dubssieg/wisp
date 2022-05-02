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
        'domain_ref': [5, 10000, 100],
        'phylum_ref': [5, 10000, 500],
        'group_ref': [5, 10000, 750],
        'order_ref': [4, 10000, 1000],
        'family_ref': [4, 10000, 1500],
        # params for for your sample here
        'domain_sample': [5, 10000, 500],
        'phylum_sample': [5, 10000, 500],
        'group_sample': [5, 10000, 750],
        'order_sample': [4, 10000, 1000],
        'family_sample': [4, 10000, 1000],
        # 'input' : location of genomes
        'input': "/udd/sidubois/Stage/Genomes/",
        # 'output' : output for database
        'output': "data/"
    }
    save_json("wisp_params", params_job)


my_params()
