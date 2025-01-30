"""Creates the XGB models"""
import os
import logging
from utils import setup_logger
from copy import copy
from xgboost import config_context, DMatrix, train
from xgboost.core import XGBoostError

logger = setup_logger(__name__, level=logging.INFO)


def make_model(output_dir: str, datas: dict, database_json: str,
               model_params: dict, taxo_level: str, taxo_target: str,) :
    """Builds the model and saves it"""
    # Creating the booster
    config_context(booster='gbtree', min_child_weight=1, tree_method='approx', predictor='cpu_predictor')
    levels = ["root", "domain", "phylum", "group", "order", "family", "specie"]
    level_up: int = levels.index(taxo_level) + 1
    next_level: str = levels[level_up]

    if not taxo_level in datas['mappings'] and not taxo_level == 'root':
        logger.error(f"Database does not contain {taxo_level} level.")
        raise ValueError

    mappings: dict = copy(datas['mappings'][next_level])
    try:
        number_taxa: int = mappings.pop('number_taxa')
    except KeyError:
        logger.error("Key 'number_taxa was missing, for unknown reason")
        number_taxa: int = len(mappings)         # Defining default value

    model_params['num_class']= number_taxa
    num_rounds_boosting = model_params.pop('num_rounds_boosting', 10)  # Récupère et supprime la clé

    temp_dir = f"{output_dir}/tmp"
    database_name = os.path.splitext(os.path.basename(database_json))[0]
    model_dir = f"{output_dir}/model/{database_name}"
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    # We will be creating temporary LibSVM files in order to make our model learn on those,
    # then destroy files in order to save space
    temp_dataset = f"{temp_dir}/{taxo_target}_{taxo_level}.txt"
    with open(temp_dataset, 'w', encoding='utf-8') as libsvm_writer:
        for data_by_genome in datas['datas']:
            if taxo_level == 'root' or data_by_genome[taxo_level] == taxo_target: # aps de root dans data_by_genome
                # Sample should be kept for model
                for read in data_by_genome['datas']:
                    id_label_next_level = mappings[data_by_genome[next_level]]
                    label_next_level = data_by_genome[next_level]
                    kmer_pairs = ' '.join([str(k) + ':' + str(v) for k, v in read.items()]) # Each read is a dict with code:count for kmer
                    line = f"{id_label_next_level} {kmer_pairs} #{label_next_level}\n" # sert pour l'eval et non xgboost
                    libsvm_writer.write(line)

    model_output_path = f"{model_dir}/{taxo_target}_{taxo_level}.json"
    config_output_path = f"{model_dir}/{taxo_target}_{taxo_level}_params.json"

    try:
        # Creating the model
        bst = train(model_params, DMatrix(temp_dataset+"?format=libsvm"), num_rounds_boosting) # Booster
        bst.save_model(model_output_path)  # Saving the model and its params        # Must go to model_dir

    except XGBoostError as e:
        logger.error("Error XgBoost ", e)    # Invalid dataset, we don't want to keep current level
        return None, None

    with open(config_output_path, 'w', encoding='utf-8') as jwriter:
        jwriter.write(bst.save_config())

    os.remove(temp_dataset)     # Destroying the temporary directory and its contents
    return model_output_path, config_output_path  # Returning the target file


